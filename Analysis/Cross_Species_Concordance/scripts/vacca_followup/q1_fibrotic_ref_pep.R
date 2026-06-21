#!/usr/bin/env Rscript
# ============================================================================
# q1_fibrotic_ref_pep.R   (Vacca follow-up Q1)
#
# Build the FIBROTIC human reference from Vacca MOESM6 (NES sheet) as the
# per-pathway mean of the two "Severe vs Mild" human columns:
#     fibrotic_ref = mean( 'UCAM/VCU: Severe vs Mild', 'EPoS: Severe vs Mild' )
# Then PEP-correlate each of our 4 mouse diets' KEGG-pathway NES profile
# (fgsea_mouse_results.csv) against that fibrotic reference over the SHARED
# KEGG-pathway set.  Sanity target (script 10): FPC ~0.92, MCD ~0.69,
# HFD ~0.64, CDAHFD ~0.64 (Pearson PEP).
#
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(BASE, "results/vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# KEGG name -> MOESM6 pathway-name space (UPPER, underscores->spaces, drop KEGG_)
pnrm <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))

cat(strrep("=", 78), "\n")
cat("Q1  Fibrotic human reference (Severe vs Mild) -- diet PEP correlation\n")
cat(strrep("=", 78), "\n")

# -- (1) Vacca human FIBROTIC reference (MOESM6 NES sheet) --------------------
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"),
                               sheet = "NES"))
setnames(p6, 1, "pathway")
p6[, pname := pnrm(pathway)]
fib_cols <- c("UCAM/VCU: Severe vs Mild", "EPoS: Severe vs Mild")
stopifnot(all(fib_cols %in% names(p6)))
fib_mat <- sapply(fib_cols, function(cc) suppressWarnings(as.numeric(p6[[cc]])))
fib_ref <- rowMeans(fib_mat, na.rm = TRUE)
n_both  <- sum(rowSums(is.finite(fib_mat)) == 2L)
human_ref <- data.table(pname = p6$pname, h_nes = fib_ref)
human_ref <- human_ref[is.finite(h_nes)]
cat(sprintf("Fibrotic ref = mean(%s).\n  MOESM6 pathways: %d total; %d finite in fibrotic ref (%d have BOTH human cohorts).\n",
            paste(fib_cols, collapse = " , "), nrow(p6), nrow(human_ref), n_both))

# -- (2) Our per-diet KEGG NES (fgsea_mouse_results.csv) ----------------------
mr <- fread(file.path(CS, "fgsea_mouse_results.csv"))[grepl("^KEGG_", pathway)]
mr[, pname := pnrm(pathway)]
diets  <- c("MCD", "HFD", "CDAHFD", "FPC")
shared <- intersect(unique(mr$pname), human_ref$pname)
cat(sprintf("Shared KEGG pathways (ours intersect Vacca fibrotic ref): %d.\n",
            length(shared)))
href_v <- human_ref[match(shared, pname)]$h_nes

# -- (3) PEP correlation per diet over shared pathways ------------------------
res <- rbindlist(lapply(diets, function(d) {
  v <- mr[source == d][match(shared, pname)]
  ok <- is.finite(v$NES) & is.finite(href_v)
  data.table(
    our_diet     = d,
    n_pathways   = sum(ok),
    pep_pearson  = suppressWarnings(cor(v$NES, href_v, method = "pearson",  use = "complete.obs")),
    pep_spearman = suppressWarnings(cor(v$NES, href_v, method = "spearman", use = "complete.obs")))
}))
setorder(res, -pep_pearson)
res[, rank_pearson := .I]

cat("\n=== Our 4 diets -- PEP correlation vs Vacca FIBROTIC reference ===\n")
print(res[, .(our_diet, n_pathways,
              pep_pearson  = round(pep_pearson, 3),
              pep_spearman = round(pep_spearman, 3),
              rank_pearson)])

cat("\nSanity vs script 10 (Pearson PEP expected): FPC ~0.92 | MCD ~0.69 | HFD ~0.64 | CDAHFD ~0.64\n")
sanity <- c(FPC = 0.92, MCD = 0.69, HFD = 0.64, CDAHFD = 0.64)
chk <- res[, .(our_diet, pep_pearson = round(pep_pearson, 3),
               expected = sanity[our_diet],
               abs_diff = round(abs(pep_pearson - sanity[our_diet]), 3))]
print(chk[order(match(our_diet, c("FPC","MCD","HFD","CDAHFD")))])
cat(sprintf("Max |observed - expected| = %.3f  -> %s\n",
            max(chk$abs_diff), if (max(chk$abs_diff) < 0.02) "REPRODUCES script 10" else "DIVERGES"))

out <- copy(res)
out[, c("expected_pearson") := sanity[our_diet]]
fwrite(out, file.path(OUT, "q1_fibrotic_ref_pep.csv"))
cat("\nWrote:", file.path(OUT, "q1_fibrotic_ref_pep.csv"), "\n")
cat(sprintf("SHARED_PATHWAY_COUNT=%d\n", length(shared)))
cat(strrep("=", 78), "\n")
