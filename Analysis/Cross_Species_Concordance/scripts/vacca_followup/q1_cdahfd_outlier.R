#!/usr/bin/env Rscript
# ============================================================================
# q1_cdahfd_outlier.R  (Vacca follow-up Q1)
#
# QUESTION: Why does our CDAHFD diet rank LOW (pathway-PEP 0.64, tied with HFD)
# while Vacca's own CDAHFD models score high? Hypothesis: our CDAHFD data
# (GSE162876) is an OUTLIER vs Vacca's three 6J-CDAHFD-F45 models.
#
# TEST: directly correlate (Pearson + Spearman, over shared KEGG pathways) our
# CDAHFD pathway-NES profile against EACH of Vacca's 3 CDAHFD models
# (6J-CDAHFD-F45-6W/8W/12W). Then quantify how different our CDAHFD is from
# theirs, contrasted against:
#   (a) within-Vacca CDAHFD model-model agreement (the "expected" similarity),
#   (b) our FPC vs Vacca's WD models (positive control: FPC is our best diet),
#   (c) our HFD vs Vacca's HFD model + cross-diet baselines (negative control).
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark")
set.seed(42)
pnrm <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))  # KEGG -> MOESM6 space
pr   <- function(a, b) suppressWarnings(cor(a, b, method = "pearson",  use = "complete.obs"))
sp   <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))

# ── Vacca per-model pathway NES (MOESM6) ─────────────────────────────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway"); p6[, pname := pnrm(pathway)]
model_cols <- setdiff(names(p6), c("pathway", "pname",
                                   grep("UCAM|EPoS", names(p6), value = TRUE)))
for (cc in model_cols) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))

# Vacca model groupings of interest
vacca_cdahfd <- c("6J-CDAHFD-F45-6W", "6J-CDAHFD-F45-8W", "6J-CDAHFD-F45-12W")
vacca_hfd    <- c("6J-HFD-F45-28W")
# WD models = the pure-Western diet group (FPC analog). Use the high-quality
# canonical WD set (exclude FG/AMLD/CCL4/STZ-combination and ALIOS/GAN variants
# for the *primary* FPC control; report the broad WD set too).
vacca_wd_core <- c("6J-WD-C0.2-32W", "MC4R-WD-C0.2-21W")
vacca_wd_all  <- grep("-WD-", model_cols, value = TRUE)   # any Western-diet model

# ── Our per-diet KEGG pathway NES ────────────────────────────────────────────
mr <- fread(file.path(CS, "fgsea_mouse_results.csv"))[grepl("^KEGG_", pathway)]
mr[, pname := pnrm(pathway)]
ourNES <- function(d) { v <- mr[source == d]; setNames(v$NES, v$pname) }

# shared KEGG pathway universe (our diets ∩ Vacca pathway table)
shared <- intersect(unique(mr$pname), p6[is.finite(rowSums(sapply(model_cols, function(c) p6[[c]]), na.rm=TRUE))]$pname)
shared <- intersect(unique(mr$pname), p6$pname)
cat(sprintf("Shared KEGG pathways (ours ∩ Vacca MOESM6): %d\n\n", length(shared)))
href_idx <- match(shared, p6$pname)

# helper: correlate OUR diet d against a Vacca model column mc over shared paths
corr_one <- function(d, mc) {
  ov <- ourNES(d)[shared]
  vv <- p6[[mc]][href_idx]
  data.table(our_diet = d, vacca_model = mc, pearson = pr(ov, vv), spearman = sp(ov, vv),
             n = sum(is.finite(ov) & is.finite(vv)))
}

# ── (1) OUR CDAHFD vs EACH Vacca CDAHFD model ────────────────────────────────
cd <- rbindlist(lapply(vacca_cdahfd, function(mc) corr_one("CDAHFD", mc)))
cat("=== OUR CDAHFD (GSE162876) vs EACH Vacca CDAHFD model (shared KEGG NES) ===\n")
print(cd)
cat(sprintf("  mean Pearson = %.3f   mean Spearman = %.3f\n\n",
            mean(cd$pearson), mean(cd$spearman)))

# ── (2) WITHIN-Vacca CDAHFD model-model agreement (the expected similarity) ───
pairs <- combn(vacca_cdahfd, 2, simplify = FALSE)
vv_cd <- rbindlist(lapply(pairs, function(pp) {
  a <- p6[[pp[1]]][href_idx]; b <- p6[[pp[2]]][href_idx]
  data.table(model_a = pp[1], model_b = pp[2], pearson = pr(a, b), spearman = sp(a, b))
}))
cat("=== WITHIN-Vacca CDAHFD model-model agreement (expected ceiling) ===\n")
print(vv_cd)
cat(sprintf("  mean within-Vacca-CDAHFD Pearson = %.3f   Spearman = %.3f\n\n",
            mean(vv_cd$pearson), mean(vv_cd$spearman)))

# ── (3) POSITIVE CONTROL: OUR FPC vs Vacca WD models ─────────────────────────
fpc_core <- rbindlist(lapply(vacca_wd_core, function(mc) corr_one("FPC", mc)))
fpc_all  <- rbindlist(lapply(vacca_wd_all,  function(mc) corr_one("FPC", mc)))
cat("=== POS CONTROL: OUR FPC vs Vacca WD-core models ===\n"); print(fpc_core)
cat(sprintf("  FPC vs WD-core  : mean Pearson = %.3f   Spearman = %.3f\n",
            mean(fpc_core$pearson), mean(fpc_core$spearman)))
cat(sprintf("  FPC vs WD-all(%d): mean Pearson = %.3f   Spearman = %.3f\n\n",
            length(vacca_wd_all), mean(fpc_all$pearson), mean(fpc_all$spearman)))

# ── (4) NEG-CONTROL baselines: our CDAHFD vs OTHER Vacca diets ────────────────
#    (is our CDAHFD closer to Vacca CDAHFD than to a random Vacca model?)
cd_vs_all <- rbindlist(lapply(model_cols, function(mc) corr_one("CDAHFD", mc)))
cd_vs_all[, is_vacca_cdahfd := vacca_model %in% vacca_cdahfd]
cat("=== Where does our CDAHFD rank across ALL Vacca models (by Pearson) ===\n")
setorder(cd_vs_all, -pearson)
cd_vs_all[, rank := .I]
print(cd_vs_all[, .(rank, vacca_model, pearson = round(pearson, 3),
                    spearman = round(spearman, 3), is_vacca_cdahfd)])
best_cd_rank <- min(cd_vs_all[is_vacca_cdahfd == TRUE]$rank)
cat(sprintf("\n  Best-ranked Vacca-CDAHFD model among ALL models for our CDAHFD: rank %d / %d\n",
            best_cd_rank, nrow(cd_vs_all)))
cat(sprintf("  mean Pearson our-CDAHFD vs {Vacca CDAHFD} = %.3f  vs {all other Vacca} = %.3f\n\n",
            mean(cd_vs_all[is_vacca_cdahfd == TRUE]$pearson),
            mean(cd_vs_all[is_vacca_cdahfd == FALSE]$pearson)))

# Same ranking diagnostic for FPC (does FPC sit ON its WD analogs?)
fpc_vs_all <- rbindlist(lapply(model_cols, function(mc) corr_one("FPC", mc)))
fpc_vs_all[, is_vacca_wd := vacca_model %in% vacca_wd_all]
setorder(fpc_vs_all, -pearson); fpc_vs_all[, rank := .I]
best_fpc_rank <- min(fpc_vs_all[is_vacca_wd == TRUE]$rank)
cat(sprintf("  (contrast) Best-ranked Vacca-WD model among ALL for our FPC: rank %d / %d\n",
            best_fpc_rank, nrow(fpc_vs_all)))

# ── Summary table ────────────────────────────────────────────────────────────
summ <- rbindlist(list(
  data.table(comparison = "our_CDAHFD_vs_Vacca_CDAHFD",   mean_pearson = mean(cd$pearson),       mean_spearman = mean(cd$spearman),       n_models = nrow(cd)),
  data.table(comparison = "within_Vacca_CDAHFD",          mean_pearson = mean(vv_cd$pearson),    mean_spearman = mean(vv_cd$spearman),    n_models = nrow(vv_cd)),
  data.table(comparison = "our_FPC_vs_Vacca_WDcore",      mean_pearson = mean(fpc_core$pearson), mean_spearman = mean(fpc_core$spearman), n_models = nrow(fpc_core)),
  data.table(comparison = "our_FPC_vs_Vacca_WDall",       mean_pearson = mean(fpc_all$pearson),  mean_spearman = mean(fpc_all$spearman),  n_models = nrow(fpc_all)),
  data.table(comparison = "our_CDAHFD_vs_other_Vacca",    mean_pearson = mean(cd_vs_all[is_vacca_cdahfd==FALSE]$pearson), mean_spearman = mean(cd_vs_all[is_vacca_cdahfd==FALSE]$spearman), n_models = sum(cd_vs_all$is_vacca_cdahfd==FALSE))
))
summ[, `:=`(mean_pearson = round(mean_pearson, 3), mean_spearman = round(mean_spearman, 3))]
cat("\n=== SUMMARY ===\n"); print(summ)

# write outputs
fwrite(cd_vs_all, file.path(OUT, "q1_cdahfd_outlier.csv"))
fwrite(summ,      file.path(OUT, "q1_cdahfd_outlier_summary.csv"))

cat("\n", strrep("=", 78), "\n")
cat(sprintf("VERDICT: our CDAHFD-Vacca-CDAHFD Pearson %.3f vs within-Vacca-CDAHFD %.3f.\n",
            mean(cd$pearson), mean(vv_cd$pearson)))
cat(sprintf("  Gap (within - ours) = %.3f. FPC-vs-WD pos-control Pearson = %.3f.\n",
            mean(vv_cd$pearson) - mean(cd$pearson), mean(fpc_core$pearson)))
cat(strrep("=", 78), "\n")
