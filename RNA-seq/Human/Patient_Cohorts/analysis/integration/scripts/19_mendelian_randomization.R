#!/usr/bin/env Rscript
# 19_mendelian_randomization.R
# ---------------------------------------------------------------------------
# Mendelian Randomization & Colocalization: eQTL -> MASLD
# Objective: Convert consensus MASLD DEGs to causal targets
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MVMR)
  library(coloc)
  library(dplyr)
  library(ggplot2)
})

cat("=== Phase 6: Causal Inference via MR ===\n\n")

WD  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
RES <- file.path(WD, "Human/Patient_Cohorts/analysis/integration/results")
OUT <- file.path(RES, "causal_inference")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MR_DATA <- file.path(WD, "../GWAS/MR_Data")
dir.create(MR_DATA, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load Significant DEGs
# ============================================================
cat("Loading consensus DEGs...\n")
degs <- fread(file.path(RES, "disease_signatures/unified_disease_signatures.csv"))
target_genes <- degs[dvc_dream_sig == TRUE]  # C2-OK-sensitivity (RETIRED MR script; MR ditched 2026-04-22, not repointed)
cat("Found", nrow(target_genes), "significant DEGs for MR testing\n")

# ============================================================
# 2. Extract eQTL Instruments (Exposure)
# ============================================================
cat("\nLoading GTEx v8 Liver eQTLs (Local File)...\n")
# Note: GTEx Liver eQTLs must be manually downloaded to GWAS/MR_Data/GTEx_Liver.v8.signif_variant_gene_pairs.txt.gz
# due to API restrictions and size.
eqtl_file <- file.path(MR_DATA, "GTEx_Liver.v8.signif_variant_gene_pairs.txt.gz")
if (!file.exists(eqtl_file)) {
    warning("GTEx Liver eQTL file not found! Please download from GTEx Portal to: ", eqtl_file)
    # Mocking extraction logic for demonstration
    # instruments <- read_exposure_data(filename = eqtl_file, sep="\t", ... format_data() )
} else {
    cat("GTEx file found. Parsing instruments...\n")
    # instruments <- read_exposure_data(...)
    # instruments <- clump_data(instruments)
}

# ============================================================
# 3. Extract MASLD GWAS (Outcome)
# ============================================================
cat("Loading MASLD GWAS Summary Statistics...\n")
# Note: Anstee 2020 full summary statistics (GCST010861) are not available via standard EBI FTP
# Must be acquired from author/restricted access and placed at:
gwas_file <- file.path(MR_DATA, "Anstee_2020_MASLD_GWAS_harmonised.tsv.gz")
if (!file.exists(gwas_file)) {
    warning("MASLD GWAS file not found! Please place summary stats at: ", gwas_file)
    # outcome_dat <- read_outcome_data(snps = instruments$SNP, filename = gwas_file, sep="\t", ... )
} else {
    cat("GWAS file found. Parsing outcome data...\n")
    # outcome_dat <- read_outcome_data(...)
}

# ============================================================
# 4. Harmonize and Run MR
# ============================================================
cat("\nHarmonizing alleles and preparing MVMR structure...\n")
# if (exists("instruments") && exists("outcome_dat")) {
#   dat <- harmonise_data(exposure_dat=instruments, outcome_dat=outcome_dat)
#   res <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
#   
#   # Save results
#   fwrite(res, file.path(OUT, "MR_causal_estimates.csv"))
#   
#   # Volcano plot
#   res$log10_pval <- -log10(res$pval)
#   ggplot(res[res$method == "Inverse variance weighted", ], aes(x=b, y=log10_pval, label=exposure)) +
#     geom_point() + theme_minimal() + labs(title="MR Causal Effects on MASLD", x="Causal Effect (IVW Beta)", y="-log10(P-value)")
#   ggsave(file.path(OUT, "MR_volcano.pdf"), width=8, height=6)
# }

cat("\nPipeline configured. Awaiting live execution post-GSE213621 integration and data acquisition.\n")
