#!/usr/bin/env Rscript
# ============================================================================
# q1_adv_fpc_vs_individual_models.R   (Q1 ADVERSARIAL)
#
# CLAIM UNDER TEST: "Western/FPC is most human-proximal" (FPC rank-1 of OUR 4
# diets at fibrotic-PEP 0.92).
#
# Adversarial re-examination:
#  (A) POOLED RANKING. Vacca's group-level PEP collapses 33 heterogeneous mouse
#      models into a handful of diet groups. Re-compute the SAME fibrotic-PEP
#      for our 4 diets and DROP THEM INTO the per-model ranking of Vacca's 33
#      individual models (q1_twoarm_pep_per_model.csv: fibro_pep_P column,
#      built against the identical fibrotic reference). Where does FPC land
#      among 33 + 4 = 37 individual mouse profiles?
#  (B) PER-GROUP. Compare FPC's fibrotic-PEP against the DISTRIBUTION of
#      individual Vacca WD/GAN/AMLN/HFD/CDHFD models -> is FPC inside or above
#      the best individual WD/GAN model?
#  (C) PATHWAY LEVERAGE. Is FPC's 0.92 driven by a few pathways?
#       - leave-one-pathway-out (LOPO) jackknife: max drop in r when 1 pathway
#         removed; how many pathways must be removed to push r below the #2 diet
#         (MCD ~0.69) and below the median Vacca model.
#       - top per-pathway contribution (z_model * z_human / (n-1)).
#       - compare LOPO fragility of FPC vs MCD/HFD/CDAHFD (is FPC MORE fragile?).
#
# Env: rnaseq.  Reference = mean(UCAM/VCU Sev-vs-Mild, EPoS Sev-vs-Mild), the
# fibrotic reference used by q1_fibrotic_ref_pep.R / script 10.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
RES   <- file.path(BASE, "results/vacca_benchmark")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

pnrm <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))
pear <- function(a, b) suppressWarnings(cor(a, b, method = "pearson",  use = "complete.obs"))

cat(strrep("=", 78), "\n")
cat("Q1 ADVERSARIAL: FPC vs INDIVIDUAL Vacca mouse models (fibrotic-PEP)\n")
cat(strrep("=", 78), "\n")

# ── (0) Human FIBROTIC reference (MOESM6 NES) ────────────────────────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway")
p6[, pname := pnrm(pathway)]
fib_cols <- c("UCAM/VCU: Severe vs Mild", "EPoS: Severe vs Mild")
stopifnot(all(fib_cols %in% names(p6)))
fib_mat <- sapply(fib_cols, function(cc) suppressWarnings(as.numeric(p6[[cc]])))
href <- data.table(pname = p6$pname, h_nes = rowMeans(fib_mat, na.rm = TRUE))
href <- href[is.finite(h_nes)]

# Vacca per-model NES columns (everything not human ref)
human_cols <- grep("UCAM|VCU|EPoS", names(p6), value = TRUE)
model_cols <- setdiff(names(p6)[-c(1)], c(human_cols, "pname"))
for (cc in model_cols) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))

# ── (1) Our 4 diets fibrotic-PEP over SHARED KEGG (identical to script 10) ────
mr <- fread(file.path(CS, "fgsea_mouse_results.csv"))[grepl("^KEGG_", pathway)]
mr[, pname := pnrm(pathway)]
diets  <- c("MCD","HFD","CDAHFD","FPC")
shared <- intersect(unique(mr$pname), href$pname)
href_v <- href[match(shared, pname)]$h_nes

ours <- rbindlist(lapply(diets, function(d) {
  v <- mr[source == d][match(shared, pname)]$NES
  data.table(model = paste0("OURS_", d), fibro_pep_P = pear(v, href_v),
             diet_group = "OURS", source = "ours", n_path = sum(is.finite(v) & is.finite(href_v)))
}))
cat(sprintf("Our diets over %d shared KEGG pathways:\n", length(shared)))
print(ours[order(-fibro_pep_P), .(model, fibro_pep_P = round(fibro_pep_P,3), n_path)])

# ── (2) Vacca per-model fibrotic-PEP (recompute from MOESM6, mouse only) ──────
# Recompute directly so ours & theirs use identical Pearson-on-full-pathway-set
vacca_per <- rbindlist(lapply(model_cols, function(mc) {
  mv <- p6[[mc]]
  data.table(model = mc, fibro_pep_P = pear(mv, href$h_nes),
             diet_group = NA_character_, source = "vacca", n_path = sum(is.finite(mv)))
}))
vacca_per <- vacca_per[!grepl("^R[-.]", model) & is.finite(fibro_pep_P)]  # drop rat
# attach published diet_group from the prior per-model table
pm <- fread(file.path(RES, "q1_twoarm_pep_per_model.csv"))
pm[, key := toupper(gsub("[^A-Za-z0-9]","",model))]
vacca_per[, key := toupper(gsub("[^A-Za-z0-9]","",model))]
vacca_per[pm, diet_group := i.diet_group, on = "key"]
vacca_per[, key := NULL]
cat(sprintf("\nVacca individual mouse models (recomputed fibrotic-PEP): %d\n", nrow(vacca_per)))

# ── (A) POOLED RANKING: our 4 diets + Vacca 33 individual models ─────────────
pool <- rbindlist(list(vacca_per, ours), use.names = TRUE)
setorder(pool, -fibro_pep_P)
pool[, rank := .I]
N <- nrow(pool)
fpc_rank <- pool[model == "OURS_FPC", rank]
fpc_pep  <- pool[model == "OURS_FPC", fibro_pep_P]
n_vacca_above_fpc <- sum(pool[source=="vacca", fibro_pep_P] > fpc_pep)
fpc_pctile <- 100 * (N - fpc_rank) / (N - 1)
vacca_med  <- median(vacca_per$fibro_pep_P)
vacca_max  <- max(vacca_per$fibro_pep_P)
best_vacca <- pool[source=="vacca"][which.max(fibro_pep_P)]

cat("\n--- (A) POOLED individual-model ranking (37 profiles) — TOP 15 ---\n")
print(pool[1:15, .(rank, model, source, diet_group, fibro_pep_P = round(fibro_pep_P,3))])
cat(sprintf("\n>> FPC pooled rank = %d / %d  (percentile %.0f%%, fibro-PEP=%.3f)\n",
            fpc_rank, N, fpc_pctile, fpc_pep))
cat(sprintf(">> # Vacca individual models with HIGHER fibrotic-PEP than FPC: %d / %d\n",
            n_vacca_above_fpc, nrow(vacca_per)))
cat(sprintf(">> Vacca individual-model fibrotic-PEP: median=%.3f, max=%.3f (best=%s, group=%s)\n",
            vacca_med, vacca_max, best_vacca$model, best_vacca$diet_group))
for (d in diets) {
  rr <- pool[model==paste0("OURS_",d)]
  cat(sprintf("   OURS_%s: pooled rank %d/%d, PEP=%.3f\n", d, rr$rank, N, rr$fibro_pep_P))
}

# ── (B) PER-GROUP: FPC vs individual WD/GAN/AMLN/HFD/CDHFD distribution ───────
cat("\n--- (B) Vacca individual-model fibrotic-PEP by diet_group ---\n")
grp <- vacca_per[!is.na(diet_group), .(n=.N, min=round(min(fibro_pep_P),3),
        median=round(median(fibro_pep_P),3), max=round(max(fibro_pep_P),3)),
        by = diet_group][order(-median)]
print(grp)
# Western-like groups (WD/GAN/AMLN diets in Vacca nomenclature)
west_like <- vacca_per[grepl("WD|GAN|AMLN|AMLD", diet_group)]
cat(sprintf("\nWestern-like Vacca models (WD/GAN/AMLN/AMLD groups): n=%d, PEP median=%.3f, max=%.3f\n",
            nrow(west_like), median(west_like$fibro_pep_P), max(west_like$fibro_pep_P)))
cat(sprintf(">> FPC PEP=%.3f vs best Western-like Vacca model=%.3f  -> FPC %s the best individual Western model\n",
            fpc_pep, max(west_like$fibro_pep_P),
            ifelse(fpc_pep > max(west_like$fibro_pep_P), "EXCEEDS", "is BELOW")))
cat(sprintf(">> FPC percentile within Western-like distribution: %.0f%%\n",
            100*mean(west_like$fibro_pep_P < fpc_pep)))

# ── (C) PATHWAY LEVERAGE: leave-one-pathway-out jackknife for each OUR diet ──
cat("\n--- (C) Pathway leverage / LOPO jackknife (is FPC's r a few-pathway artifact?) ---\n")
lopo_stats <- rbindlist(lapply(diets, function(d) {
  v  <- mr[source == d][match(shared, pname)]$NES
  ok <- is.finite(v) & is.finite(href_v)
  vv <- v[ok]; hh <- href_v[ok]; pn <- shared[ok]
  r_full <- cor(vv, hh)
  # LOPO: r with each single pathway removed
  r_drop <- sapply(seq_along(vv), function(i) cor(vv[-i], hh[-i]))
  # per-pathway contribution to Pearson r (standardized cross-product)
  zc <- scale(vv)[,1]; zh <- scale(hh)[,1]
  contrib <- zc * zh / (length(vv) - 1)               # sums to r
  ord <- order(-contrib)
  top5_frac <- sum(contrib[ord[1:5]]) / r_full        # frac of r from top-5 pathways
  # how many top pathways removed to push r below MCD baseline (0.69) and below vacca median
  r_after_removing_top_k <- function(k) {
    keep <- setdiff(seq_along(vv), ord[1:k]); cor(vv[keep], hh[keep])
  }
  k_below_mcd <- NA_integer_; k_below_vmed <- NA_integer_
  for (k in 1:min(40, length(vv)-3)) {
    rk <- r_after_removing_top_k(k)
    if (is.na(k_below_mcd)  && rk < 0.69)      k_below_mcd  <- k
    if (is.na(k_below_vmed) && rk < vacca_med) k_below_vmed <- k
  }
  data.table(diet = d, n_path = length(vv), r_full = r_full,
             lopo_min_r = min(r_drop), lopo_max_r = max(r_drop),
             lopo_range = max(r_drop) - min(r_drop),
             top5_contrib_frac = top5_frac,
             top_pathway = pn[ord[1]], top_contrib = contrib[ord[1]],
             k_top_to_drop_below_0.69 = k_below_mcd,
             k_top_to_drop_below_vacca_median = k_below_vmed)
}))
print(lopo_stats[, .(diet, n_path, r_full = round(r_full,3),
                     lopo_range = round(lopo_range,4),
                     top5_contrib_frac = round(top5_contrib_frac,3),
                     k_top_to_drop_below_0.69, k_top_to_drop_below_vacca_median)])

cat("\nTop-5 leverage pathways for FPC:\n")
{ d <- "FPC"
  v  <- mr[source == d][match(shared, pname)]$NES
  ok <- is.finite(v) & is.finite(href_v)
  vv <- v[ok]; hh <- href_v[ok]; pn <- shared[ok]
  zc <- scale(vv)[,1]; zh <- scale(hh)[,1]
  contrib <- zc*zh/(length(vv)-1); ord <- order(-contrib)
  print(data.table(pathway = pn[ord[1:5]],
                   contrib = round(contrib[ord[1:5]],4),
                   fpc_NES = round(vv[ord[1:5]],2),
                   human_NES = round(hh[ord[1:5]],2)))
}

# ── verdict numbers ──────────────────────────────────────────────────────────
fpc_lopo <- lopo_stats[diet=="FPC"]
robust <- (n_vacca_above_fpc <= 3) &&            # FPC near the very top of individuals
          (fpc_lopo$top5_contrib_frac < 0.5) &&  # not a top-5-pathway artifact
          (is.na(fpc_lopo$k_top_to_drop_below_0.69) || fpc_lopo$k_top_to_drop_below_0.69 >= 5)

cat("\n", strrep("=", 78), "\n", sep="")
cat(sprintf("VERDICT INPUTS:\n  FPC pooled rank %d/%d (%.0f pctile); %d Vacca individuals beat FPC\n",
            fpc_rank, N, fpc_pctile, n_vacca_above_fpc))
cat(sprintf("  FPC top-5-pathway frac of r = %.2f ; LOPO range = %.4f ; need to drop top-%s pathways to fall below MCD(0.69)\n",
            fpc_lopo$top5_contrib_frac, fpc_lopo$lopo_range,
            ifelse(is.na(fpc_lopo$k_top_to_drop_below_0.69),">40",as.character(fpc_lopo$k_top_to_drop_below_0.69))))
cat(sprintf("  => 'FPC most human-proximal' is %s\n", ifelse(robust, "ROBUST", "FRAGILE/ARTIFACT")))
cat(strrep("=", 78), "\n")

# ── write outputs ────────────────────────────────────────────────────────────
fwrite(pool[, .(rank, model, source, diet_group, fibro_pep_P, n_path)],
       file.path(RES, "q1_adv_pooled_individual_ranking.csv"))
fwrite(lopo_stats, file.path(RES, "q1_adv_fpc_pathway_leverage.csv"))
summ <- data.table(
  metric = c("fpc_fibro_pep","fpc_pooled_rank","pool_n","fpc_percentile",
             "n_vacca_individuals_above_fpc","vacca_individual_median","vacca_individual_max",
             "best_vacca_model","best_vacca_group","west_like_n","west_like_max_pep",
             "fpc_exceeds_best_western","fpc_top5_contrib_frac","fpc_lopo_range",
             "k_top_drop_below_mcd","k_top_drop_below_vacca_median","verdict_robust"),
  value  = c(round(fpc_pep,3), fpc_rank, N, round(fpc_pctile,1),
             n_vacca_above_fpc, round(vacca_med,3), round(vacca_max,3),
             best_vacca$model, best_vacca$diet_group, nrow(west_like), round(max(west_like$fibro_pep_P),3),
             fpc_pep > max(west_like$fibro_pep_P), round(fpc_lopo$top5_contrib_frac,3),
             round(fpc_lopo$lopo_range,4),
             ifelse(is.na(fpc_lopo$k_top_to_drop_below_0.69),NA,fpc_lopo$k_top_to_drop_below_0.69),
             ifelse(is.na(fpc_lopo$k_top_to_drop_below_vacca_median),NA,fpc_lopo$k_top_to_drop_below_vacca_median),
             robust))
fwrite(summ, file.path(RES, "q1_adv_fpc_vs_individual_models.csv"))
cat("Wrote:\n  ", file.path(RES,"q1_adv_fpc_vs_individual_models.csv"), "\n  ",
    file.path(RES,"q1_adv_pooled_individual_ranking.csv"), "\n  ",
    file.path(RES,"q1_adv_fpc_pathway_leverage.csv"), "\n")
