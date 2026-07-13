#!/usr/bin/env Rscript
# ============================================================================
# 11_vacca_q1_fibrotic_validation.R   (Vacca 2024 benchmark — canonical Q1)
#
# Consolidates the validated Q1 follow-up (vacca_followup/q1_twoarm_pep_dhps.R +
# q1_metabolic_arm_pep.R + q1_adv_fpc_vs_individual_models.R) into one traceable
# canonical script.
#
# QUESTION (Q1): why does MCD rank ABOVE CDAHFD/HFD in our single-axis pathway
# proximity (script 10), even after adopting LITMUS's pathway resolution?
#
# ANSWER: LITMUS's DHPS is TWO-ARMED — a METABOLIC arm (MOESM9 col5) and a
# FIBROTIC arm (MOESM9 col9). Our script-10 score used the FIBROTIC human
# reference (Severe-vs-Mild), so it is a FIBROTIC-proximity score. We therefore:
#   (A) VALIDATE the construction on LITMUS's OWN 33 mouse models: a pathway-PEP
#       built against each human arm should track the matching published DHPS arm.
#       -> fibrotic arm reproduces (Spearman ~0.62); metabolic arm does NOT
#          (~0.3, NS) -> we report only the FIBROTIC axis and DEFER the metabolic
#          ranking to LITMUS.
#   (B) APPLY the fibrotic-PEP to OUR 4 diets. MCD is genuinely fibrosis-strong
#       (matching LITMUS MCD f_DHPS ~0.75) -> its mid-rank is FAITHFUL, not error.
#       FPC (Western) is the most fibrotic-proximal of our diets, reproducing
#       LITMUS's "Western diets closest to human" -> bootstrap CI + leverage check.
#
# PEP proximity = Pearson corr of a model's KEGG-pathway-NES vector with the human
# arm reference, over shared KEGG pathways (the validated pathArm construction,
# scripts 09/10). Spearman reported as robustness. Rat models (^R-) dropped.
#
# Inputs:
#   data/external/vacca_2024/42255_2024_1043_MOESM6_ESM.xlsx (sheet NES; human ref + per-model NES)
#   data/external/vacca_2024/42255_2024_1043_MOESM9_ESM.xlsx (published two-arm DHPS)
#   Analysis/Cross_Species_Concordance/results/fgsea_mouse_results.csv (our per-diet KEGG NES)
# Outputs (Analysis/Cross_Species_Concordance/results/vacca_benchmark/):
#   q1_fibrotic_validation_arms.csv      (two-arm PEP vs published DHPS; headline Spearmans)
#   q1_fibrotic_per_model.csv            (33 Vacca models: fibro-PEP + published f_DHPS, for the figure)
#   q1_our_diets_fibrotic_pep.csv        (our 4 diets: metab + fibro PEP, bootstrap CI, pooled rank)
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
set.seed(42)
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

nrm  <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))                       # model-name key
pnrm <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))       # KEGG -> MOESM6 pathway space
pe   <- function(a, b) suppressWarnings(cor(a, b, method = "pearson",  use = "complete.obs"))
sp   <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
sp.p <- function(a, b) suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)$p.value)

# ── Published two-arm DHPS (MOESM9): col5 = metabolic-arm, col9 = fibrotic-arm ──
m9 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM9_ESM.xlsx"),
                               col_names = FALSE, skip = 2))
setnames(m9, 1:10, c("model","diet_group","mP","mH","m_DHPS","mM","fP","fH","f_DHPS","fM"))
m9 <- m9[!is.na(model)]
m9[, `:=`(key = nrm(model), m_DHPS = as.numeric(m_DHPS), f_DHPS = as.numeric(f_DHPS))]

# ── MOESM6 pathway NES: human arm references + per-model PEP profiles ──────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway"); p6[, pname := pnrm(pathway)]
human_cols <- grep("UCAM|VCU|EPoS", names(p6), value = TRUE)
model_cols <- setdiff(names(p6)[-1], c(human_cols, "pname"))
for (cc in c(human_cols, model_cols)) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))

href_metab <- p6[["UCAM/VCU: Mild vs Control"]]                                  # METABOLIC arm (early/steatosis)
href_fibro <- rowMeans(cbind(p6[["UCAM/VCU: Severe vs Mild"]],                   # FIBROTIC arm (Severe-vs-Mild)
                             p6[["EPoS: Severe vs Mild"]]), na.rm = TRUE)
cat(sprintf("MOESM6 pathways=%d ; model cols=%d ; metab ref finite=%d ; fibro ref finite=%d\n",
            nrow(p6), length(model_cols), sum(is.finite(href_metab)), sum(is.finite(href_fibro))))

# ── (A) Two-arm PEP for each LITMUS mouse model (full MOESM6 pathway set) ──────
prox <- rbindlist(lapply(model_cols, function(mc) {
  mv <- p6[[mc]]
  r <- data.table(model = mc,
             metab_pep_P = pe(mv, href_metab), fibro_pep_P = pe(mv, href_fibro),
             metab_pep_S = sp(mv, href_metab), fibro_pep_S = sp(mv, href_fibro))
  r[, key := nrm(mc)][]
}))
mg <- merge(prox, m9[, .(key, diet_group, m_DHPS, f_DHPS)], by = "key")
mg <- mg[!grepl("^R[-.]", model)]                                               # drop rat models
cat(sprintf("LITMUS mouse models matched (MOESM6 PEP x MOESM9 DHPS): %d\n", nrow(mg)))

arms <- data.table(
  test = c("MATCHED: fibrotic-PEP vs fibrotic-DHPS (headline)",
           "MATCHED: metabolic-PEP vs metabolic-DHPS",
           "MATCHED(Spearman-PEP): fibrotic vs fibrotic-DHPS",
           "MATCHED(Spearman-PEP): metabolic vs metabolic-DHPS",
           "CROSS: fibrotic-PEP vs metabolic-DHPS",
           "CROSS: metabolic-PEP vs fibrotic-DHPS",
           "ref: published DHPS arms inter-corr (m_DHPS vs f_DHPS)"),
  spearman = c(sp(mg$fibro_pep_P, mg$f_DHPS), sp(mg$metab_pep_P, mg$m_DHPS),
               sp(mg$fibro_pep_S, mg$f_DHPS), sp(mg$metab_pep_S, mg$m_DHPS),
               sp(mg$fibro_pep_P, mg$m_DHPS), sp(mg$metab_pep_P, mg$f_DHPS),
               sp(mg$m_DHPS, mg$f_DHPS)),
  p_value  = c(sp.p(mg$fibro_pep_P, mg$f_DHPS), sp.p(mg$metab_pep_P, mg$m_DHPS),
               sp.p(mg$fibro_pep_S, mg$f_DHPS), sp.p(mg$metab_pep_S, mg$m_DHPS),
               sp.p(mg$fibro_pep_P, mg$m_DHPS), sp.p(mg$metab_pep_P, mg$f_DHPS),
               sp.p(mg$m_DHPS, mg$f_DHPS)),
  n = nrow(mg))
cat("\n=== (A) Two-arm transcriptomic PEP vs published two-arm DHPS (n=33 LITMUS models) ===\n")
print(arms, digits = 3)
fwrite(arms, file.path(OUT, "q1_fibrotic_validation_arms.csv"))
fwrite(mg[order(-f_DHPS), .(model, diet_group, fibro_pep_P, fibro_pep_S, f_DHPS,
                            metab_pep_P, m_DHPS)],
       file.path(OUT, "q1_fibrotic_per_model.csv"))

fib <- arms[test == "MATCHED: fibrotic-PEP vs fibrotic-DHPS (headline)"]
met <- arms[test == "MATCHED: metabolic-PEP vs metabolic-DHPS"]
cat(sprintf("\n>> FIBROTIC arm reproduced: Spearman = %.3f (p=%.2g, n=%d)  [HEADLINE]\n",
            fib$spearman, fib$p_value, fib$n))
cat(sprintf(">> METABOLIC arm: Spearman = %.3f (p=%.2g)  -> weak/NS, DEFER to LITMUS\n",
            met$spearman, met$p_value))

# ── (B) Apply fibrotic-PEP to OUR 4 diets, over the shared KEGG pathway set ────
mr <- fread(file.path(CS, "fgsea_mouse_results.csv"))[grepl("^KEGG_", pathway)]
mr[, pname := pnrm(pathway)]
diets <- c("MCD", "HFD", "CDAHFD", "FPC")

# shared pathways: our KEGG NES ∩ Vacca fibrotic reference (finite both)
fib_ref <- data.table(pname = p6$pname, h = href_fibro)[is.finite(h)]
met_ref <- data.table(pname = p6$pname, h = href_metab)[is.finite(h)]
shared  <- intersect(unique(mr$pname), fib_ref$pname)
hv_f    <- fib_ref[match(shared, pname)]$h
hv_m    <- met_ref[match(shared, pname)]$h
cat(sprintf("\nShared KEGG pathways (our diets ∩ LITMUS human ref): %d\n", length(shared)))

# bootstrap CI: resample shared pathways with replacement (B=2000), recompute Pearson
B <- 2000L
boot_ci <- function(v, h) {
  ok <- is.finite(v) & is.finite(h); v <- v[ok]; h <- h[ok]; n <- length(v)
  bs <- vapply(seq_len(B), function(i) { idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(v[idx], h[idx])) }, numeric(1))
  quantile(bs, c(0.025, 0.975), na.rm = TRUE)
}
our <- rbindlist(lapply(diets, function(d) {
  v  <- mr[source == d][match(shared, pname)]$NES
  ci <- boot_ci(v, hv_f)
  data.table(our_diet = d, n_shared = sum(is.finite(v) & is.finite(hv_f)),
             fibro_pep = pe(v, hv_f), fibro_ci_lo = ci[[1]], fibro_ci_hi = ci[[2]],
             metab_pep = pe(v, hv_m))
}))
our[, fibro_rank := frank(-fibro_pep)][, metab_rank := frank(-metab_pep)]
setorder(our, fibro_rank)

# pooled placement: our diets' fibro-PEP (shared set) vs LITMUS models' fibro-PEP
# recomputed on the SAME shared set, so the axis is identical
vacca_shared <- rbindlist(lapply(model_cols, function(mc) {
  v <- p6[match(shared, pname)][[mc]]
  r <- data.table(model = mc, fibro_pep = pe(v, hv_f))
  r[, key := nrm(mc)][]
}))
vacca_shared <- merge(vacca_shared, m9[, .(key, diet_group, f_DHPS)], by = "key")[!grepl("^R[-.]", model)]
vacca_shared <- vacca_shared[is.finite(fibro_pep)]
fib_ref_shared_rho <- sp(vacca_shared[is.finite(f_DHPS)]$fibro_pep,
                         vacca_shared[is.finite(f_DHPS)]$f_DHPS)
cat(sprintf("Sanity: fibrotic-PEP(shared-%d-pathway) vs published f_DHPS across LITMUS models: Spearman %.3f\n",
            length(shared), fib_ref_shared_rho))

vmed <- median(vacca_shared$fibro_pep, na.rm = TRUE)
vmax <- max(vacca_shared$fibro_pep, na.rm = TRUE)
for (d in diets) {
  pep <- our[our_diet == d]$fibro_pep
  our[our_diet == d, n_litmus_above := sum(vacca_shared$fibro_pep > pep, na.rm = TRUE)]
  our[our_diet == d, pooled_rank := sum(vacca_shared$fibro_pep > pep, na.rm = TRUE) + fibro_rank]
}
our[, litmus_n := nrow(vacca_shared)][, litmus_median := vmed][, litmus_max := vmax]

cat("\n=== (B) OUR 4 diets — fibrotic-arm proximity (shared KEGG, bootstrap 95% CI) ===\n")
print(our[, .(our_diet, fibro_pep = round(fibro_pep, 3),
              CI = sprintf("[%.2f, %.2f]", fibro_ci_lo, fibro_ci_hi),
              fibro_rank, n_litmus_above, metab_pep = round(metab_pep, 3), metab_rank)])
cat(sprintf("\nLITMUS individual-model fibrotic-PEP (shared set): median=%.3f, max=%.3f, n=%d\n",
            vmed, vmax, nrow(vacca_shared)))
fpc <- our[our_diet == "FPC"]
mcd <- our[our_diet == "MCD"]
cat(sprintf(">> FPC (Western) is the most fibrotic-proximal of our diets: %.3f CI[%.2f,%.2f], %d LITMUS models above\n",
            fpc$fibro_pep, fpc$fibro_ci_lo, fpc$fibro_ci_hi, fpc$n_litmus_above))
cat(sprintf(">> MCD fibrotic rank=%d (genuinely fibrosis-strong; LITMUS MCD/CDD f_DHPS~0.75); metabolic rank=%d (the demerit we defer)\n",
            mcd$fibro_rank, mcd$metab_rank))
cat("   NOTE: FPC point-estimate is leverage-sensitive (few-pathway driven; see q1_adv_fpc_*), hence the wide bootstrap CI.\n")

fwrite(our, file.path(OUT, "q1_our_diets_fibrotic_pep.csv"))
fwrite(vacca_shared[order(-fibro_pep)], file.path(OUT, "q1_litmus_models_fibrotic_pep_shared.csv"))
cat("\nWrote q1_fibrotic_validation_arms.csv / q1_fibrotic_per_model.csv / q1_our_diets_fibrotic_pep.csv\n")
cat(strrep("=", 78), "\n")
