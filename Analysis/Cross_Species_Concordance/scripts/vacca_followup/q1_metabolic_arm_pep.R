#!/usr/bin/env Rscript
# ============================================================================
# q1_metabolic_arm_pep.R
#
# Q1 (Vacca follow-up): Build the METABOLIC human pathway reference =
#   MOESM6 'UCAM/VCU: Mild vs Control' KEGG pathway-NES (steatosis/early axis,
#   i.e. mild MASLD vs healthy control). PEP-correlate (Pearson AND Spearman
#   over SHARED KEGG pathways) each of OUR 4 mouse diets' KEGG pathway-NES
#   profiles (fgsea_mouse_results.csv, source==diet, names normalized to the
#   MOESM6 pathway-name space) against this metabolic reference.
#
# QUESTION: does MCD score LOW on the metabolic-arm reference (it should --
#   MCD is a lean/choline-deficient fibrotic model, NOT a metabolic-obesity
#   model), reversing the MCD>CDAHFD/HFD ordering we saw against the FIBROTIC
#   (Severe-vs-Mild) reference in script 10?
#
# For contrast we also recompute the FIBROTIC reference (Severe vs Mild,
#   UCAM/VCU + EPoS mean, exactly as script 10) so the metabolic-vs-fibrotic
#   reordering of our diets is read off side by side.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

pnrm <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))  # KEGG->MOESM6 space
pe   <- function(a, b) suppressWarnings(cor(a, b, method = "pearson",  use = "complete.obs"))
sp   <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))

# ── Vacca human pathway NES (MOESM6) ─────────────────────────────────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway"); p6[, pname := pnrm(pathway)]

num <- function(col) suppressWarnings(as.numeric(p6[[col]]))

# METABOLIC arm = Mild vs Control (steatosis/early). UCAM/VCU is the cohort with
# a healthy-Control arm (EPoS has no 'vs Control' contrast -> only UCAM/VCU here).
metab_ref  <- data.table(pname = p6$pname, h_nes = num("UCAM/VCU: Mild vs Control"))
# FIBROTIC arm = Severe vs Mild (mean of UCAM/VCU + EPoS), identical to script 10.
fib_ref    <- data.table(pname = p6$pname,
  h_nes = rowMeans(cbind(num("UCAM/VCU: Severe vs Mild"), num("EPoS: Severe vs Mild")), na.rm = TRUE))

# ── OUR per-diet KEGG pathway NES ────────────────────────────────────────────
mr <- fread(file.path(CS, "fgsea_mouse_results.csv"))[grepl("^KEGG_", pathway)]
mr[, pname := pnrm(pathway)]
diets <- c("MCD", "HFD", "CDAHFD", "FPC")

pep_against <- function(ref, label) {
  shared <- intersect(unique(mr$pname), ref[is.finite(h_nes)]$pname)
  hv <- ref[match(shared, pname)]$h_nes
  res <- rbindlist(lapply(diets, function(d) {
    v <- mr[source == d][match(shared, pname)]$NES
    data.table(arm = label, our_diet = d,
               n_shared = sum(is.finite(v) & is.finite(hv)),
               pearson  = pe(v, hv),
               spearman = sp(v, hv))
  }))
  cat(sprintf("\n[%s arm] shared KEGG pathways (ours x Vacca ref) = %d\n", label, length(shared)))
  res
}

metab <- pep_against(metab_ref, "metabolic_MildVsControl")
fibro <- pep_against(fib_ref,    "fibrotic_SevereVsMild")

cat("\n=== METABOLIC arm (Mild vs Control) -- our 4 diets, ranked by Pearson ===\n")
print(metab[order(-pearson)], digits = 3)
cat("\n  ranked by Spearman:\n")
print(metab[order(-spearman), .(our_diet, spearman = round(spearman, 3))])

cat("\n=== FIBROTIC arm (Severe vs Mild) -- our 4 diets, ranked by Pearson (script-10 reference) ===\n")
print(fibro[order(-pearson)], digits = 3)

# side-by-side: metabolic vs fibrotic proximity per diet + rank shift
cmp <- merge(
  metab[, .(our_diet, metab_pearson = pearson, metab_spearman = spearman)],
  fibro[, .(our_diet, fibro_pearson = pearson, fibro_spearman = spearman)],
  by = "our_diet")
cmp[, metab_rank := frank(-metab_pearson)]
cmp[, fibro_rank := frank(-fibro_pearson)]
setorder(cmp, metab_rank)

cat("\n=== SIDE-BY-SIDE: metabolic-arm vs fibrotic-arm proximity (Pearson) ===\n")
print(cmp, digits = 3)

mcd_metab_rank <- cmp[our_diet == "MCD"]$metab_rank
mcd_fibro_rank <- cmp[our_diet == "MCD"]$fibro_rank
cat(sprintf("\nMCD metabolic-arm rank = %d/4 (Pearson %.3f, Spearman %.3f)\n",
            mcd_metab_rank, cmp[our_diet=="MCD"]$metab_pearson, cmp[our_diet=="MCD"]$metab_spearman))
cat(sprintf("MCD fibrotic-arm rank  = %d/4 (Pearson %.3f) -- script-10 axis\n",
            mcd_fibro_rank, cmp[our_diet=="MCD"]$fibro_pearson))
cat(sprintf("=> MCD %s on the metabolic reference relative to fibrotic.\n",
            ifelse(mcd_metab_rank > mcd_fibro_rank, "DROPS (scores LOWER, as expected for a lean model)",
                   ifelse(mcd_metab_rank == mcd_fibro_rank, "is UNCHANGED", "RISES (unexpected)"))))

# write tidy long output
allres <- rbind(metab, fibro)
fwrite(allres, file.path(OUT, "q1_metabolic_arm_pep.csv"))
fwrite(cmp,    file.path(OUT, "q1_metabolic_vs_fibrotic_pep.csv"))
cat("\nWrote:", file.path(OUT, "q1_metabolic_arm_pep.csv"), "\n")
cat("Wrote:", file.path(OUT, "q1_metabolic_vs_fibrotic_pep.csv"), "\n")
