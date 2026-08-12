# Shared data loader for MASLD publication figures
# Centralized Ensembl->symbol mapping + column normalization
# Source AFTER publication_theme.R

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INTEGRATION <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
INT_RESULTS <- file.path(INTEGRATION, "results/integration")
INT_META    <- file.path(INTEGRATION, "metadata")
CAUSAL      <- file.path(BASE, "RNA-seq/results/causal_inference")
DRUG        <- file.path(BASE, "RNA-seq/results/drug_repurposing")
GWAS_RES    <- file.path(BASE, "RNA-seq/results/gwas_spatial_convergence")
ME          <- file.path(BASE, "RNA-seq/results/multi_evidence")
CONCORDANCE <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
MOUSE       <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results")
VALIDATION  <- file.path(BASE, "RNA-seq/results/validation")
AUDIT       <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
PER_STUDY   <- file.path(INTEGRATION, "results/per_study")
SIGS        <- file.path(INTEGRATION, "results/disease_signatures")
FIG_OUT     <- Sys.getenv("MASLD_FIGURE_OUTPUT_ROOT", file.path(BASE, "figures"))
PROTEOMICS_DIR <- file.path(BASE, "Analysis/Proteomics/results")
ATAC_DIR       <- file.path(BASE, "Analysis/ATAC/Human_Multiome")
SPATIAL_DIR    <- file.path(BASE, "Analysis/Spatial/results")
PROGRAM_CONTEXT_DIR <- file.path(BASE, "Analysis/Multimodal_Program_Projection/results")
NETWORK_DIR    <- file.path(BASE, "RNA-seq/results/network")

# ═══════════════════════════════════════════════════════════════════════════════
# GWAS registry ancestry / trait map — SINGLE SOURCE OF TRUTH (added 2026-07-05)
# The 50-GWAS COLOC portfolio (23 legacy + 27 MVP strata; MVP integrated 2026-07-04)
# enumerates ancestry + trait per study in gwas_registry.tsv. Fig 2 panels MUST
# derive ancestry via gwas_ancestry() — NOT the legacy grepl() heuristic
# (BBJ->EAS / PanUKBB_AFR->AFR / PanUKBB_CSA->SAS / else->EUR), which silently
# misroutes EVERY MVP stratum (MVP_*_AMR/AFR/EAS/EUR) into EUR and has no AMR bin.
# Verified 2026-07-05: all 50 susie_coloc_all_gwas `gwas_name` values match a
# registry `study_name` exactly, so the strict join below never spuriously errors.
# ═══════════════════════════════════════════════════════════════════════════════
GWAS_REGISTRY <- file.path(BASE, "GWAS/finemapping/config/gwas_registry.tsv")
# Canonical ancestry factor order + colourblind-safe palette (EUR reference first,
# then alphabetical). AMR added 2026-07-05 with MVP. Python panels that cannot
# source this file should mirror these hex values.
GWAS_ANCESTRY_LEVELS <- c("EUR", "AFR", "AMR", "EAS", "SAS")
ANCESTRY_COLORS <- c(EUR = "#4C72B0", AFR = "#55A868", AMR = "#DD8452",
                     EAS = "#C44E52", SAS = "#8172B3")

.gwas_registry_cache <- NULL
load_gwas_registry <- function() {
  if (!is.null(.gwas_registry_cache)) return(.gwas_registry_cache)
  if (!file.exists(GWAS_REGISTRY)) { message("WARNING: ", GWAS_REGISTRY, " not found"); return(NULL) }
  .gwas_registry_cache <<- fread(GWAS_REGISTRY)
  .gwas_registry_cache
}

# study_name -> ancestry via the registry. Strict by default: any unmapped name is
# an ERROR (a new cohort must be registered, never silently defaulted to EUR).
gwas_ancestry <- function(gwas_name, strict = TRUE) {
  reg <- load_gwas_registry()
  m <- reg$ancestry[match(as.character(gwas_name), reg$study_name)]
  if (any(is.na(m))) {
    miss <- unique(as.character(gwas_name)[is.na(m)])
    msg <- sprintf("gwas_ancestry(): %d GWAS name(s) not in registry: %s",
                   length(miss), paste(head(miss, 20), collapse = ", "))
    if (strict) stop(msg) else message("WARNING: ", msg)
  }
  factor(m, levels = GWAS_ANCESTRY_LEVELS)
}

# study_name -> trait token (registry trait_type is only binary/quantitative, so the
# trait label is parsed from the name). Covers legacy + MVP naming; longest/most
# specific tokens first so e.g. "ALT" cannot pre-empt a longer match. Returns NA
# (not a silent "NAFLD") for unmatched names.
gwas_trait <- function(gwas_name) {
  g <- as.character(gwas_name)
  out <- rep(NA_character_, length(g))
  toks <- c("ChronLiver", "Cirrhosis", "CHIRHEP", "Albumin", "Platelet", "PDFF",
            "NAFLD", "NASH", "HCC", "ALT", "AST", "GGT", "cirrhosis")
  for (t in toks) {
    hit <- is.na(out) & grepl(t, g, fixed = TRUE)
    out[hit] <- t
  }
  # normalize synonyms to display labels
  out[out %in% c("ChronLiver")]            <- "Chronic liver disease"
  out[out %in% c("CHIRHEP", "cirrhosis")]  <- "Cirrhosis"
  out
}

# ═══════════════════════════════════════════════════════════════════════════════
# Cohort metadata edge-case guards — SINGLE SOURCE OF TRUTH (added 2026-06-24)
# Any per-sample fibrosis_stage / NAS binning MUST route through staged_disease_meta()
# / clean_nas_meta() (defined below) so contaminated cohorts can't silently enter
# stage-stratified plots. The DEG/GSEA per-stage results (load_fibrosis_stage_dream
# etc., from Script 14b) are already clean — these guards are for figures/analyses
# that read RAW unified_metadata.csv and bin it.
#   - GSE213621 (Chen): COARSE fibrosis — the study reports grouped bins (F0F1/F2/F3F4)
#     and 00_harmonize_metadata.R maps F0F1->1, F2->2, F3F4->3. NOT true Kleiner F0-F4
#     (no true F0/F4); largest cohort (~1/3 of all F1/F2/F3 disease samples).
#   - PRJNA512027 (Gerhard): dropped 2026-05-15 (L0/S0 library-prep batch confound).
#   - GSE135251/213621/240729: sex is INFERRED (XIST/DDX3Y k-means), not annotated.
# STAGED_FIB_COHORTS mirrors the 14b DE allowlist (the 6 true-Kleiner cohorts).
# ═══════════════════════════════════════════════════════════════════════════════
COARSE_STAGE_COHORTS <- c("GSE213621")
DROPPED_COHORTS      <- c("PRJNA512027")
INFERRED_SEX_COHORTS <- c("GSE135251", "GSE213621", "GSE240729")
STAGED_FIB_COHORTS   <- c("GSE130970", "GSE135251", "GSE162694",
                          "GSE174478", "GSE193066", "GSE240729")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure output directories — SINGLE SOURCE OF TRUTH
# All figure scripts MUST use these constants. NEVER hardcode output paths.
# To add a new directory, update this file AND figures/README.md.
# ═══════════════════════════════════════════════════════════════════════════════
FIG_MAIN  <- file.path(FIG_OUT, "main")
FIG_SUPP  <- file.path(FIG_OUT, "supplementary")
FIG_MISC  <- file.path(FIG_OUT, "misc")
FIG_ARCHIVE <- file.path(FIG_OUT, "archive")

# Deprecated deconvolution-attribution diagnostics retained for provenance.
# These are not manuscript-facing figure directories.
FIG_DECONV_C2_ARCHIVE   <- file.path(FIG_ARCHIVE, "c2_deconv_raw_reference")
FIG_DECONV_RECT_ARCHIVE <- file.path(FIG_ARCHIVE, "rectangle_deconv_sensitivity")

# Consolidated methods-validation / robustness / QC parent (2026-06-04).
# Ten formerly-scattered validation/robustness/QC supp dirs now live as
# subfolders here. The 5 repointed constants below (FIGS01_DIR qc_validation,
# FIGS_SENS_DIR sensitivity, FIGS_HCAUDIT_DIR healthy_control_audit,
# FIGS_BATCH_DIR batch_correction, FIGS_QUANT_DIR quantification) keep their
# names so consumer scripts need no edit; the 5 below are new.
FIGS_METHVAL_DIR   <- file.path(FIG_SUPP, "figS_methods_validation")
FIGS_MEGAVAL_DIR   <- file.path(FIGS_METHVAL_DIR, "mega_validation")
FIGS_MULTIMETH_DIR <- file.path(FIGS_METHVAL_DIR, "multimethod_validation")
FIGS_INTVAL_DIR    <- file.path(FIGS_METHVAL_DIR, "integration_value")
FIGS_ROBUST_DIR    <- file.path(FIGS_METHVAL_DIR, "robustness")
FIGS_LFCSENS_DIR   <- file.path(FIGS_METHVAL_DIR, "lfc_sensitivity")

# Main-figure directory compatibility constants. Final semantic placement is
# governed by docs/PAPER.md and docs/ROADMAP.md, not these legacy names.
FIG1_DIR  <- file.path(FIG_MAIN, "fig1_atlas_overview")           # Atlas + cohorts
FIG2_DIR  <- file.path(FIG_MAIN, "fig3_RNAseq")                   # Final Fig 3: cross-sectional established-state transcriptomics; FIG2_DIR retained only for back-compat.
FIG3_DIR  <- file.path(FIG_MAIN, "fig2_genetics")                 # Genetics / GWAS-eQTL (main Fig 2; dir renamed fig3_regulatory_architecture -> fig2_genetics 2026-06-12; constant name FIG3_DIR kept for back-compat across ~21 consumer scripts). NB: this is main Fig 2, distinct from FIG2_DIR (fig3_RNAseq) above.
FIG4_DIR  <- file.path(FIG_MAIN, "fig4_validation")               # Final Fig 4: assay-native molecular and physical context; directory name is compatibility-only.
FIG5_DIR  <- file.path(FIG_MAIN, "fig5_convergence")              # Final Fig 5: MASLD Gene Catalog; directory name is compatibility-only.

# Supplementary figures (S1-S10 + sensitivity + therapeutics)
FIGS01_DIR    <- file.path(FIGS_METHVAL_DIR, "qc_validation")  # was figS01_qc_validation (consolidated 2026-06-04)
FIGS02_DIR    <- file.path(FIG_SUPP, "figS02_progression")
FIGS03_DIR    <- file.path(FIG_SUPP, "figS03_deconvolution")
FIGS04_DIR    <- file.path(FIG_SUPP, "figS04_coloc")
FIGS05_DIR    <- file.path(FIG_SUPP, "figS05_epigenomic_spatial")
FIGS06_DIR    <- file.path(FIG_SUPP, "figS06_cross_species")
FIGS07_DIR    <- file.path(FIG_SUPP, "figS07_ncrna")
FIGS08_DIR    <- file.path(FIG_SUPP, "figS08_subtyping_convergence")
FIGS09_DIR    <- file.path(FIG_SUPP, "figS09_multi_ancestry")
FIGS10_DIR    <- file.path(FIG_SUPP, "figS10_prediction")
FIGS_SENS_DIR <- file.path(FIGS_METHVAL_DIR, "sensitivity")  # was figS_sensitivity (consolidated 2026-06-04)
FIGS_HEPSUB_DIR <- file.path(FIG_SUPP, "figS_hepatocyte_subtypes")
FIGS_THERA_DIR  <- file.path(FIG_SUPP, "figS_therapeutics")  # demoted fig6 panels
FIGS_NET_DIR    <- file.path(FIG_SUPP, "figS_network")       # Bayesian multiplex network
FIGS_RORA_DIR   <- file.path(FIG_SUPP, "figS_rora_case_study")  # RORA multi-modal case study
FIGS_CELLTYPE_DIR <- file.path(FIG_SUPP, "figS_celltype_biology")  # Cell-type-resolved MASLD biology (27 analyses; A-L themes)
FIGS_MCP_DIR      <- file.path(FIG_SUPP, "figS_mcp")            # Multi-cellular programs (cNMF + DIALOGUE)
FIGS_CONV_EVID_DIR <- file.path(FIG_SUPP, "figS_convergence_evidence")  # Convergence evidence score (Script 46d)
FIGS_HCAUDIT_DIR  <- file.path(FIGS_METHVAL_DIR, "healthy_control_audit")  # Healthy-control audit; was figS_healthy_control_audit (consolidated 2026-06-04)
FIGS_SEX_DIR      <- file.path(FIG_SUPP, "figS_sex_dimorphism")  # Sex-dimorphic biology demoted from Fig 2 (2026-04-29)
FIGS_BATCH_DIR    <- file.path(FIGS_METHVAL_DIR, "batch_correction")  # Harmony batch-correction adequacy; was figS_batch_correction (consolidated 2026-06-04)
FIGS_GRANULAR_DIR <- file.path(FIG_SUPP, "figS_granular_staging")  # Two-transition decomposition + F3 sub-state multi-modal (2026-05-06)
FIGS_SCDRS_DIR    <- file.path(FIG_SUPP, "figS_scdrs")           # Consolidated scDRS supp figs (bulk-DEG anchor + GWAS-anchored; 2026-05-12)
FIGS_SCDRS_DATA_DIR <- file.path(FIGS_SCDRS_DIR, "panel_data")   # Per-panel CSVs for caption transparency
FIGS_STAGECCC_DIR <- file.path(FIG_SUPP, "stage_ccc")            # Stage-stratified CCC trajectories (Scripts 349/349b/349c/351/352/353; 2026-05-13)
FIGS_PROTEO_DIR   <- file.path(FIG_SUPP, "figS_proteomics")       # Proteomics supp (DIA-MS volcano/enrichment + liver->blood decoupling demoted from Fig 4, 2026-07-02)
FIGS_GLP1RA_DIR   <- file.path(FIG_SUPP, "figS_glp1ra")          # GLP-1RA / incretin axis (mechanism-only; demoted from Fig 5, 2026-07-13)

# Hotspot autocorrelation modules (Pipeline 14; 2026-05-17). FLAT layout —
# all composites land directly under figS_hotspot/; all individual panels under
# figS_hotspot/panels/ with topic prefixes. legacy/ holds retired 509 outputs.
FIGS_HOTSPOT_DIR         <- file.path(FIG_SUPP, "figS_hotspot")
FIGS_HOTSPOT_PANELS_DIR  <- file.path(FIGS_HOTSPOT_DIR, "panels")
FIGS_HOTSPOT_DATA_DIR    <- file.path(FIGS_HOTSPOT_PANELS_DIR, "data")

# Separate future Cas13-screen-paper figures. This constant is not part of the
# standalone Resource manuscript or its release graph.
FIGS_CAS13LIB_DIR <- file.path(BASE, "Cas13_Library_Design", "figures")

# Quantification comparison (STAR vs Kallisto sensitivity)
FIGS_QUANT_DIR <- file.path(FIGS_METHVAL_DIR, "quantification")  # was figS_quantification (consolidated 2026-06-04)

# Create all directories
for (d in c(FIG_MAIN, FIG_SUPP, FIG_MISC, FIG_ARCHIVE,
            FIG_DECONV_C2_ARCHIVE, FIG_DECONV_RECT_ARCHIVE,
            FIG1_DIR, FIG2_DIR, FIG3_DIR, FIG4_DIR, FIG5_DIR,
            FIGS01_DIR, FIGS02_DIR, FIGS03_DIR, FIGS04_DIR, FIGS05_DIR,
            FIGS06_DIR, FIGS07_DIR, FIGS08_DIR, FIGS09_DIR, FIGS10_DIR,
            FIGS_SENS_DIR, FIGS_HEPSUB_DIR, FIGS_THERA_DIR,
            FIGS_NET_DIR, FIGS_RORA_DIR, FIGS_CELLTYPE_DIR, FIGS_MCP_DIR,
            FIGS_CONV_EVID_DIR, FIGS_HCAUDIT_DIR, FIGS_SEX_DIR, FIGS_BATCH_DIR,
            FIGS_GRANULAR_DIR, FIGS_SCDRS_DIR, FIGS_SCDRS_DATA_DIR,
            FIGS_STAGECCC_DIR, FIGS_PROTEO_DIR, FIGS_GLP1RA_DIR,
            FIGS_HOTSPOT_DIR, FIGS_HOTSPOT_PANELS_DIR, FIGS_HOTSPOT_DATA_DIR,
            FIGS_CAS13LIB_DIR,
            FIGS_QUANT_DIR,
            FIGS_METHVAL_DIR, FIGS_MEGAVAL_DIR, FIGS_MULTIMETH_DIR,
            FIGS_INTVAL_DIR, FIGS_ROBUST_DIR, FIGS_LFCSENS_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---------------------------------------------------------------------------
# Gene map: Ensembl -> HGNC symbol
# ---------------------------------------------------------------------------
.gene_map_cache <- NULL

load_gene_map <- function() {
  if (!is.null(.gene_map_cache)) return(.gene_map_cache)
  atlas_path <- file.path(ME, "multi_evidence_atlas.csv")
  if (!file.exists(atlas_path)) {
    message("WARNING: multi_evidence_atlas.csv not found; gene map empty")
    .gene_map_cache <<- data.table(ensembl_clean = character(), symbol = character())
    return(.gene_map_cache)
  }
  atlas <- fread(atlas_path, select = c("ensembl_id", "human_symbol"))
  atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]
  atlas <- atlas[!is.na(human_symbol) & human_symbol != "", .(ensembl_clean, symbol = human_symbol)]
  atlas <- unique(atlas, by = "ensembl_clean")

  # Supplement from concordance atlas
  conc_path <- file.path(CONCORDANCE, "concordance_atlas_unified.csv")
  if (file.exists(conc_path)) {
    conc <- fread(conc_path, select = c("mouse_gene_id", "human_symbol"))
    conc <- conc[!is.na(human_symbol) & human_symbol != ""]
    conc[, ensembl_clean := sub("\\..*", "", mouse_gene_id)]
    conc <- conc[, .(ensembl_clean, symbol = human_symbol)]
    conc <- conc[!ensembl_clean %in% atlas$ensembl_clean]
    atlas <- rbind(atlas, unique(conc, by = "ensembl_clean"))
  }
  .gene_map_cache <<- atlas
  atlas
}

# Add symbol column to data.table with Ensembl gene IDs
add_symbols <- function(dt, gene_col = "gene") {
  if ("symbol" %in% names(dt)) return(dt)   # already present
  gm <- load_gene_map()
  dt[, ensembl_clean := sub("\\..*", "", get(gene_col))]
  dt <- merge(dt, gm, by = "ensembl_clean", all.x = TRUE)
  dt[is.na(symbol), symbol := ensembl_clean]
  dt[, ensembl_clean := NULL]
  dt
}

# ---------------------------------------------------------------------------
# Data loaders (cached)
# ---------------------------------------------------------------------------

# --- Merged DGEList (for per-patient analyses) ---
.dge_cache <- NULL
load_merged_dge <- function() {
  if (!is.null(.dge_cache)) return(.dge_cache)
  f <- file.path(INT_RESULTS, "merged_dge.rds")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dge <- readRDS(f)
  .dge_cache <<- dge
  dge
}

# --- Dream results ---
.dream_cache <- NULL
load_dream_results <- function() {
  if (!is.null(.dream_cache)) return(.dream_cache)
  # Canonical DEG table (hard cutover 2026-06-08): canonical_deg_results.csv.
  # Schema: gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr, AveExpr, symbol.
  # C2 sweep (2026-06-08): output column names are emitted as bulk_* (matching
  # the multi-evidence atlas S1 columns) so canonical-DEG and atlas figures share
  # one naming. The retired dream_* aliases were dropped.
  f <- file.path(INT_RESULTS, "canonical_deg_results.csv")
  if (!file.exists(f)) {
    message("WARNING: ", f, " not found")
    return(NULL)
  }
  dt <- fread(f)
  # Map ashr columns to the bulk_* names
  if ("shrunk_logFC" %in% names(dt) && !"bulk_shrunk_logFC" %in% names(dt))
    setnames(dt, "shrunk_logFC", "bulk_shrunk_logFC")
  if ("lfsr" %in% names(dt) && !"bulk_lfsr" %in% names(dt))
    setnames(dt, "lfsr", "bulk_lfsr")
  # Normalize canonical unprefixed column names to bulk_*
  if ("padj" %in% names(dt) && !"bulk_padj" %in% names(dt))
    setnames(dt, "padj", "bulk_padj")
  if ("logFC" %in% names(dt) && !"bulk_logFC" %in% names(dt))
    setnames(dt, "logFC", "bulk_logFC")
  dt <- add_symbols(dt, "gene")
  .dream_cache <<- dt
  dt
}

# --- DEG classification helper ---
# Returns a logical vector: TRUE if gene passes the canonical DEG threshold.
# Canonical (2026-06-29): TREAT FDR < 0.05 at lfc = 0.25 (real limma::treat, computed in
# 05h). treat() tests H0:|true logFC|<=lfc, so the effect-size floor is folded INTO the
# test -- there is NO separate |logFC|/|shrunk| filter. Supersedes the ashr lfsr+|shrunk|>0.3
# gate (those columns are retained as reference). `dt` must carry a treat_fdr column;
# accepts bulk_treat_fdr (atlas / consensus) or treat_fdr (canonical loader output).
is_dream_deg <- function(dt) {
  tf_col <- intersect(c("bulk_treat_fdr", "treat_fdr"), names(dt))[1]
  if (is.na(tf_col))
    stop("is_dream_deg: need treat_fdr (bulk_treat_fdr/treat_fdr); have: ",
         paste(names(dt), collapse = ", "))
  !is.na(dt[[tf_col]]) & dt[[tf_col]] < 0.05
}

# --- MASH vs MASL results ---
.mash_masl_cache <- NULL
load_mash_vs_masl_results <- function() {
  if (!is.null(.mash_masl_cache)) return(.mash_masl_cache)
  f <- file.path(SIGS, "nafl_vs_nash_dream.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  # Normalize to unprefixed canonical names (MASH-vs-MASL contrast; not the main
  # disease-vs-control bulk channel). Source became limma-voom under the C2 swap
  # (2026-06-08): it ships `padj` (was adj.P.Val) + `logFC`/`shrunk_logFC`/`lfsr`.
  if (!"padj" %in% names(dt)) {
    padj_src <- intersect(c("adj.P.Val"), names(dt))[1]
    if (!is.na(padj_src)) setnames(dt, padj_src, "padj")
  }
  if (!"logFC" %in% names(dt)) {
    lfc_src <- intersect(c("shrunk_logFC"), names(dt))[1]
    if (!is.na(lfc_src)) setnames(dt, lfc_src, "logFC")
  }
  dt <- add_symbols(dt, "gene")
  .mash_masl_cache <<- dt
  dt
}

# --- Meta-analysis results ---
.meta_cache <- NULL
load_meta_results <- function() {
  if (!is.null(.meta_cache)) return(.meta_cache)
  f <- file.path(INT_RESULTS, "meta_analysis_results.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  dt <- add_symbols(dt, "gene")
  .meta_cache <<- dt
  dt
}

# --- Consensus DEGs (adds tier column) ---
.consensus_cache <- NULL
load_consensus_degs <- function() {
  if (!is.null(.consensus_cache)) return(.consensus_cache)
  f <- file.path(INT_RESULTS, "consensus_degs.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  dt <- add_symbols(dt, "gene")

  # If tier column already exists, use it

  if ("tier" %in% names(dt) || "consensus_tier" %in% names(dt)) {
    if ("consensus_tier" %in% names(dt)) setnames(dt, "consensus_tier", "tier")
    .consensus_cache <<- dt
    return(dt)
  }

  # Otherwise compute tier from dream + meta
  meta <- load_meta_results()
  if (!is.null(meta)) {
    meta_slim <- meta[, .(gene, meta_logFC, meta_padj)]
    dt <- merge(dt, meta_slim, by = "gene", all.x = TRUE)
  }

  # Determine significance flags (consensus_degs.csv now carries bulk_* cols)
  padj_thr <- 0.1; lfc_thr <- 0.5
  bulk_padj_col <- intersect(c("bulk_padj", "padj"), names(dt))[1]
  bulk_lfc_col  <- intersect(c("bulk_logFC", "logFC"), names(dt))[1]

  if (!is.na(bulk_padj_col) && !is.na(bulk_lfc_col)) {
    dt[, bulk_sig_tier := get(bulk_padj_col) < padj_thr & abs(get(bulk_lfc_col)) > lfc_thr]
  } else {
    dt[, bulk_sig_tier := FALSE]
  }
  if ("meta_padj" %in% names(dt) && "meta_logFC" %in% names(dt)) {
    dt[, meta_sig_tier := meta_padj < padj_thr & abs(meta_logFC) > lfc_thr]
    # Direction concordance
    if (!is.na(bulk_lfc_col)) {
      dt[, dir_concordant := sign(get(bulk_lfc_col)) == sign(meta_logFC) | is.na(meta_logFC)]
    } else {
      dt[, dir_concordant := TRUE]
    }
  } else {
    dt[, meta_sig_tier := FALSE]
    dt[, dir_concordant := TRUE]
  }

  # Assign tiers
  dt[, tier := "Not_Consensus"]
  dt[bulk_sig_tier == TRUE & meta_sig_tier == TRUE & dir_concordant == TRUE,
     tier := "Tier1_HighConfidence"]
  dt[tier == "Not_Consensus" &
     (bulk_sig_tier == TRUE | meta_sig_tier == TRUE) &
     dir_concordant == TRUE,
     tier := "Tier2_Moderate"]
  # Tier3 needs per-study counts (simplified: use bulk_sig from consensus file)
  if ("bulk_sig" %in% names(dt)) {
    # If n_studies column exists, use it
    if ("n_studies_sig" %in% names(dt)) {
      dt[tier == "Not_Consensus" & n_studies_sig >= 3, tier := "Tier3_Exploratory"]
    }
  }
  # Clean up temp columns
  dt[, c("bulk_sig_tier", "meta_sig_tier", "dir_concordant") := NULL]
  .consensus_cache <<- dt
  dt
}

# --- Concordance atlas ---
.concordance_cache <- NULL
load_concordance_atlas <- function() {
  if (!is.null(.concordance_cache)) return(.concordance_cache)
  f <- file.path(CONCORDANCE, "concordance_atlas_unified.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if (!"symbol" %in% names(dt) && "human_symbol" %in% names(dt))
    dt[, symbol := human_symbol]
  .concordance_cache <<- dt
  dt
}

# --- Deconvolution attribution ---
.deconv_cache <- NULL
load_deconv_attribution <- function() {
  if (!is.null(.deconv_cache)) return(.deconv_cache)
  f <- file.path(CAUSAL, "deconv_attribution_scores.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  dt <- add_symbols(dt, "gene")
  # Normalize category names for color matching
  dt[, category_clean := gsub(" ", "_", category)]
  .deconv_cache <<- dt
  dt
}

# --- Sex DEG classification ---
.sex_cache <- NULL
load_sex_classification <- function() {
  if (!is.null(.sex_cache)) return(.sex_cache)
  # Prefer v3/v6-shimmed csv (canonical as of 2026-05-16 v6 cross-pillar
  # consensus rebuild — 06i shim overlays v6 onto v3 schema for downstream
  # consumers; v5 shim 06h overlay is preserved under *_v5legacy aliases).
  # Fall back to legacy sex_deg_classification.csv if shim path is missing.
  f_v3 <- file.path(INT_RESULTS, "sex_v3", "sex_deg_classification_v3.csv")
  f_lg <- file.path(INT_RESULTS, "sex_deg_classification.csv")
  f <- if (file.exists(f_v3)) f_v3 else f_lg
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  # Adapter: legacy callers expect logFC_M / logFC_F (or meta_logFC_M/F).
  # v3 csv carries beta_F / beta_M from dream M2. Alias if needed.
  if (!"logFC_M" %in% names(dt) && "beta_M" %in% names(dt))
    dt[, logFC_M := beta_M]
  if (!"logFC_F" %in% names(dt) && "beta_F" %in% names(dt))
    dt[, logFC_F := beta_F]
  if (!"meta_logFC_M" %in% names(dt) && "beta_M" %in% names(dt))
    dt[, meta_logFC_M := beta_M]
  if (!"meta_logFC_F" %in% names(dt) && "beta_F" %in% names(dt))
    dt[, meta_logFC_F := beta_F]
  # v6 cross-pillar aliases (post-2026-05-16 meta-review): expose the v6
  # consensus tier + diagnostic emit columns under stable names so fig2
  # panels can overlay a robustness flag onto the v5 fixed-effect display
  # without column-name munging. Columns are populated when 06i shim has
  # run; otherwise initialized to NA so callers can detect absence.
  for (col in c("class_v6_consensus", "class_v6_combined", "consensus_quality",
                "n_pillars_used", "stability_v6", "padj_int_5k_random",
                "q_emp_per_gene", "I2_pct", "loco_rho")) {
    if (!col %in% names(dt)) dt[, (col) := NA]
  }
  dt <- add_symbols(dt, "gene")
  .sex_cache <<- dt
  dt
}

# --- Causal inference: MR ---
# RETIRED 2026-04-22 — MR ditched from paper. load_mr_results() now returns
# NULL so callers that still reference it (e.g. figS_sensitivity.R) short-circuit
# their MR panels without breaking. Delete call sites rather than the function
# stub if you encounter residual references.
load_mr_results <- function(gwas = "ghodsian") {
  message("load_mr_results(): RETIRED 2026-04-22 — MR ditched from paper. Returning NULL.")
  NULL
}

# --- Causal inference: TWAS ---
.twas_cache <- NULL
load_twas_results <- function(gwas = "ghodsian") {
  if (!is.null(.twas_cache)) return(.twas_cache)
  f <- file.path(CAUSAL, gwas, "twas_spredixcan_liver.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  # gene_name is already symbol
  if ("gene_name" %in% names(dt)) dt[, symbol := gene_name]
  .twas_cache <<- dt
  dt
}

# --- Causal inference: TWAS combined (multi-GWAS) ---
.twas_combined_cache <- NULL
load_twas_combined <- function() {
  if (!is.null(.twas_combined_cache)) return(.twas_combined_cache)
  f <- file.path(CAUSAL, "twas_multi_gwas_combined.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("gene_name" %in% names(dt)) dt[, symbol := gene_name]
  .twas_combined_cache <<- dt
  dt
}

# --- Causal inference: COLOC ---
.coloc_cache <- NULL
load_coloc_results <- function(gwas = "ghodsian") {
  if (!is.null(.coloc_cache)) return(.coloc_cache)
  f <- file.path(CAUSAL, gwas, "coloc_results.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  dt <- add_symbols(dt, "gene")
  .coloc_cache <<- dt
  dt
}

# --- LINCS annotated compounds ---
.lincs_cache <- NULL
load_lincs_compounds <- function() {
  if (!is.null(.lincs_cache)) return(.lincs_cache)
  f <- file.path(DRUG, "lincs_annotated_compounds.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  # Merge moa.x and moa.y into single moa column
  if ("moa.x" %in% names(dt) && "moa.y" %in% names(dt)) {
    dt[, moa.x := as.character(moa.x)]
    dt[, moa.y := as.character(moa.y)]
    dt[, moa := ifelse(!is.na(moa.x) & nchar(moa.x) > 0, moa.x,
                       ifelse(!is.na(moa.y) & nchar(moa.y) > 0, moa.y, NA_character_))]
    dt[, c("moa.x", "moa.y") := NULL]
  }
  if ("target.x" %in% names(dt) && "target.y" %in% names(dt)) {
    dt[, target.x := as.character(target.x)]
    dt[, target.y := as.character(target.y)]
    dt[, target := ifelse(!is.na(target.x) & nchar(target.x) > 0, target.x,
                          ifelse(!is.na(target.y) & nchar(target.y) > 0, target.y, NA_character_))]
    dt[, c("target.x", "target.y") := NULL]
  }
  .lincs_cache <<- dt
  dt
}

# --- Network proximity ---
.netprox_cache <- NULL
load_network_proximity <- function() {
  if (!is.null(.netprox_cache)) return(.netprox_cache)
  f <- file.path(DRUG, "network_proximity", "network_proximity_scores.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .netprox_cache <<- dt
  dt
}

# --- CGP reversal hits ---
.cgp_cache <- NULL
load_cgp_hits <- function() {
  if (!is.null(.cgp_cache)) return(.cgp_cache)
  f <- file.path(DRUG, "cgp_reversal_hits.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .cgp_cache <<- dt
  dt
}

# --- Convergent drug targets ---
.drug_targets_cache <- NULL
load_drug_targets <- function() {
  if (!is.null(.drug_targets_cache)) return(.drug_targets_cache)
  # Try convergent first, then multi_layer
  f <- file.path(DRUG, "convergent_drug_targets.csv")
  if (!file.exists(f)) f <- file.path(DRUG, "multi_layer_drug_targets.csv")
  if (!file.exists(f)) { message("WARNING: drug targets file not found"); return(NULL) }
  dt <- fread(f)
  .drug_targets_cache <<- dt
  dt
}

# --- ComBat-seq sensitivity ---
# Repointed 2026-07-04 (round-2 audit D2j) from Script 34 (dream, mis-scoped over
# the full merged set) to Script 34c2 (`combat_seq_c2/`): C2-canonical LVQW engine
# on the SAME 5 control-bearing cohorts (n=846) the canonical DE was fit on
# (r=0.901 vs canonical; batch removal does not strip disease signal). Function
# names kept for back-compat; they now serve the 34c2 (LVQW C2) outputs.
.combat_dream_cache <- NULL
load_combat_dream <- function() {
  if (!is.null(.combat_dream_cache)) return(.combat_dream_cache)
  f <- file.path(AUDIT, "combat_seq_c2", "c2_combat_seq_results.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  # Normalize adj.P.Val -> padj
  if ("adj.P.Val" %in% names(dt) && !"padj" %in% names(dt))
    setnames(dt, "adj.P.Val", "padj")
  dt <- add_symbols(dt, "gene")
  .combat_dream_cache <<- dt
  dt
}

.combat_concordance_cache <- NULL
load_combat_concordance <- function() {
  if (!is.null(.combat_concordance_cache)) return(.combat_concordance_cache)
  f <- file.path(AUDIT, "combat_seq_c2", "combat_seq_concordance_c2.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  # 34c2 uses C2-canonical metric names; remap to the schema the figure expects
  # (sig_primary/sig_combat/sig_overlap/jaccard). direction_concordance and
  # lfc_pearson_r already match.
  if ("metric" %in% names(dt)) {
    remap <- c(sig_canon_p01 = "sig_primary", sig_combat_p01 = "sig_combat",
               sig_overlap_p01 = "sig_overlap", jaccard_p01 = "jaccard")
    hit <- dt$metric %in% names(remap)
    dt$metric[hit] <- remap[dt$metric[hit]]
  }
  .combat_concordance_cache <<- dt
  dt
}

# --- Variance partition ---
.vp_cache <- NULL
load_variance_partition <- function() {
  if (!is.null(.vp_cache)) return(.vp_cache)
  f <- file.path(INT_RESULTS, "variance_partition.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  dt <- add_symbols(dt, "gene")
  .vp_cache <<- dt
  dt
}

# --- Per-study DE results ---
.per_study_cache <- NULL
load_per_study_de <- function() {
  if (!is.null(.per_study_cache)) return(.per_study_cache)
  # PRJNA512027 (Gerhard 2018) removed: L0/S0 confound — excluded from project
  datasets <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE167523",
                "GSE174478", "GSE193066", "GSE213621", "GSE240729")
  all_de <- rbindlist(lapply(datasets, function(ds) {
    f <- file.path(PER_STUDY, paste0(ds, "_de_results.csv"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    if (!"dataset" %in% names(dt)) dt[, dataset := ds]
    # Normalize adj.P.Val -> padj
    if ("adj.P.Val" %in% names(dt) && !"padj" %in% names(dt))
      setnames(dt, "adj.P.Val", "padj")
    dt
  }), fill = TRUE)
  if (nrow(all_de) > 0) all_de <- add_symbols(all_de, "gene")
  .per_study_cache <<- all_de
  all_de
}

# --- Per-study summary ---
load_per_study_summary <- function() {
  f <- file.path(PER_STUDY, "per_study_summary.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  fread(f)
}

# --- Positive controls ---
.posctrl_cache <- NULL
load_positive_controls <- function() {
  if (!is.null(.posctrl_cache)) return(.posctrl_cache)
  f <- file.path(VALIDATION, "positive_control_validation.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if (!"symbol" %in% names(dt) && "gene" %in% names(dt)) {
    # gene column might already be symbols in this file
    if (!any(grepl("^ENSG", dt$gene))) {
      dt[, symbol := gene]
    } else {
      dt <- add_symbols(dt, "gene")
    }
  }
  .posctrl_cache <<- dt
  dt
}

# --- Sex-divergent permissive ---
.sex_perm_cache <- NULL
load_sex_divergent_perm <- function() {
  if (!is.null(.sex_perm_cache)) return(.sex_perm_cache)
  f <- file.path(AUDIT, "sex_divergent_permissive.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  dt <- add_symbols(dt, "gene")
  .sex_perm_cache <<- dt
  dt
}

# --- Conserved core enrichment ---
load_conserved_enr <- function() {
  f <- file.path(AUDIT, "conserved_enrichment.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  fread(f)
}

# --- CGP liver metabolic filtered ---
load_cgp_liver_filtered <- function() {
  f <- file.path(AUDIT, "cgp_liver_metabolic_filtered.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  fread(f)
}

# --- Multi-evidence atlas ---
.me_cache <- NULL
load_multi_evidence <- function() {
  if (!is.null(.me_cache)) return(.me_cache)
  f <- file.path(ME, "multi_evidence_atlas.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if (!"symbol" %in% names(dt) && "human_symbol" %in% names(dt))
    dt[, symbol := human_symbol]
  .me_cache <<- dt
  dt
}

# --- Atlas layer correlations ---
load_layer_correlations <- function() {
  f <- file.path(ME, "atlas_layer_correlations.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  fread(f)
}

# --- Benchmark results ---
load_benchmark_results <- function() {
  f <- file.path(ME, "benchmark_results.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  fread(f)
}

# --- PR curve data ---
load_pr_curve <- function() {
  f <- file.path(ME, "pr_curve_data.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  fread(f)
}

# --- GWAS convergence evidence cards ---
load_gwas_evidence <- function() {
  f <- file.path(GWAS_RES, "gwas_convergence_evidence_cards.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  fread(f)
}

# --- Unified metadata ---
load_metadata <- function() {
  f <- file.path(INT_META, "unified_metadata.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  fread(f)
}

# --- Sex-stratified meta results ---
load_sex_meta <- function() {
  f <- file.path(INT_RESULTS, "sex_stratified_meta_results.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  dt <- add_symbols(dt, "gene")
  dt
}

# --- Attribution cross-method comparison ---
load_attribution_comparison <- function() {
  f <- file.path(CAUSAL, "attribution_cross_method_comparison.csv")
  if (!file.exists(f)) {
    # Try alternative path
    f <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/deconvolution/attribution_cross_method_comparison.csv")
  }
  if (!file.exists(f)) { message("WARNING: attribution comparison not found"); return(NULL) }
  fread(f)
}

# --- Human UMAP coordinates ---
.human_umap_cache <- NULL
load_human_umap <- function() {
  if (!is.null(.human_umap_cache)) return(.human_umap_cache)
  f <- file.path(INT_RESULTS, "umap_coordinates.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .human_umap_cache <<- dt
  dt
}

# --- Mouse UMAP coordinates ---
.mouse_umap_cache <- NULL
load_mouse_umap <- function() {
  if (!is.null(.mouse_umap_cache)) return(.mouse_umap_cache)
  f <- file.path(MOUSE, "umap_coordinates.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .mouse_umap_cache <<- dt
  dt
}

# --- sc-eQTL multi-cell-type COLOC summary ---
.sceqtl_coloc_cache <- NULL
load_sceqtl_coloc <- function() {
  if (!is.null(.sceqtl_coloc_cache)) return(.sceqtl_coloc_cache)
  sceqtl_dir <- file.path(CAUSAL, "sceqtl")
  # Per-cell-type COLOC results
  cell_types <- c("hepatocyte", "endothelial_cell", "cholangiocyte", "stellate_cell")
  all_coloc <- rbindlist(lapply(cell_types, function(ct) {
    f <- file.path(sceqtl_dir, paste0("coloc_results_", ct, ".csv"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    if (!"cell_type" %in% names(dt)) dt[, cell_type := ct]
    dt
  }), fill = TRUE)
  # Also load summary
  summary_f <- file.path(sceqtl_dir, "sceqtl_multi_celltype_summary.csv")
  if (file.exists(summary_f)) {
    attr(all_coloc, "summary") <- fread(summary_f)
  }
  .sceqtl_coloc_cache <<- all_coloc
  all_coloc
}

# --- sc-eQTL MR results ---
.sceqtl_mr_cache <- NULL
load_sceqtl_mr <- function() {
  if (!is.null(.sceqtl_mr_cache)) return(.sceqtl_mr_cache)
  f <- file.path(CAUSAL, "sceqtl", "mr_results.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .sceqtl_mr_cache <<- dt
  dt
}

# --- Zenodo precomputed COLOC ---
.zenodo_coloc_cache <- NULL
load_zenodo_coloc <- function() {
  if (!is.null(.zenodo_coloc_cache)) return(.zenodo_coloc_cache)
  f <- file.path(CAUSAL, "sceqtl", "zenodo_coloc_integrated.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .zenodo_coloc_cache <<- dt
  dt
}

# --- Zenodo perturbation validation ---
.zenodo_perturb_cache <- NULL
load_zenodo_perturbation <- function() {
  if (!is.null(.zenodo_perturb_cache)) return(.zenodo_perturb_cache)
  f <- file.path(CAUSAL, "sceqtl", "zenodo_perturbation_validation.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .zenodo_perturb_cache <<- dt
  dt
}

# --- Broadaway COLOC (aggregated across all GWAS: Ghodsian, UKBB ALT, AST, GGT, PDFF) ---
.broadaway_coloc_cache <- NULL
load_broadaway_coloc <- function() {
  if (!is.null(.broadaway_coloc_cache)) return(.broadaway_coloc_cache)
  broadaway_dirs <- c("broadaway", "broadaway_ukbb", "broadaway_ukbb_ast",
                       "broadaway_ukbb_ggt", "broadaway_pdff")
  parts <- list()
  for (d in broadaway_dirs) {
    f <- file.path(CAUSAL, d, "coloc_results.csv")
    if (file.exists(f)) {
      tmp <- fread(f)
      if (nrow(tmp) > 0 && "PP.H4" %in% names(tmp)) {
        tmp <- add_symbols(tmp, "gene")
        parts[[d]] <- tmp
      }
    }
  }
  if (length(parts) == 0) { message("WARNING: No Broadaway COLOC results found"); return(NULL) }
  dt <- rbindlist(parts, fill = TRUE)
  # Keep best PP.H4 per gene across all GWAS
  dt <- dt[order(-PP.H4)][!duplicated(gene)]
  .broadaway_coloc_cache <<- dt
  dt
}

# --- ieQTL disease-interaction genes ---
.ieqtl_cache <- NULL
load_ieqtl <- function() {
  if (!is.null(.ieqtl_cache)) return(.ieqtl_cache)
  f <- file.path(CAUSAL, "sceqtl", "ieqtl_disease_genes.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .ieqtl_cache <<- dt
  dt
}

# --- UKBB ALT sc-eQTL COLOC ---
.ukbb_coloc_cache <- NULL
load_sceqtl_coloc_ukbb <- function() {
  if (!is.null(.ukbb_coloc_cache)) return(.ukbb_coloc_cache)
  sceqtl_dir <- file.path(CAUSAL, "sceqtl")
  cell_types <- c("hepatocyte", "endothelial_cell", "cholangiocyte", "stellate_cell")
  all_coloc <- rbindlist(lapply(cell_types, function(ct) {
    f <- file.path(sceqtl_dir, paste0("coloc_results_ukbb_alt_", ct, ".csv"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    if (!"cell_type" %in% names(dt)) dt[, cell_type := ct]
    dt
  }), fill = TRUE)
  summary_f <- file.path(sceqtl_dir, "sceqtl_coloc_ukbb_alt_summary.csv")
  if (file.exists(summary_f)) {
    attr(all_coloc, "summary") <- fread(summary_f)
  }
  .ukbb_coloc_cache <<- all_coloc
  all_coloc
}

# --- Fibrosis pairwise DE ---
.fib_pw_cache <- NULL
load_fibrosis_pairwise <- function() {
  if (!is.null(.fib_pw_cache)) return(.fib_pw_cache)
  f <- file.path(SIGS, "fibrosis_pairwise.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("gene" %in% names(dt)) dt <- add_symbols(dt, "gene")
  .fib_pw_cache <<- dt
  dt
}

# --- Fibrosis ordinal slopes (meta across studies) ---
.fib_slopes_cache <- NULL
load_fibrosis_slopes <- function() {
  if (!is.null(.fib_slopes_cache)) return(.fib_slopes_cache)
  f <- file.path(SIGS, "fibrosis_slopes_meta.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("gene" %in% names(dt)) dt <- add_symbols(dt, "gene")
  .fib_slopes_cache <<- dt
  dt
}

# --- Fibrosis early vs late ---
.fib_el_cache <- NULL
load_fibrosis_early_late <- function() {
  if (!is.null(.fib_el_cache)) return(.fib_el_cache)
  f <- file.path(SIGS, "fibrosis_early_late.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("gene" %in% names(dt)) dt <- add_symbols(dt, "gene")
  .fib_el_cache <<- dt
  dt
}

# --- NAS component DE ---
.nas_cache <- NULL
load_nas_components <- function() {
  if (!is.null(.nas_cache)) return(.nas_cache)
  f <- file.path(SIGS, "nas_components_all.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("gene" %in% names(dt)) dt <- add_symbols(dt, "gene")
  .nas_cache <<- dt
  dt
}

# --- Per-NAS-score dream results (Script 14b) ---
.nas_score_prog <- NULL
load_nas_score_progression <- function() {
  if (!is.null(.nas_score_prog)) return(.nas_score_prog)
  f <- file.path(SIGS, "nas_score_dream.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("adj.P.Val" %in% names(dt) && !"padj" %in% names(dt))
    setnames(dt, "adj.P.Val", "padj")
  dt <- add_symbols(dt, "gene")
  .nas_score_prog <<- dt
  dt
}

# --- Per-NAS-score sample sizes ---
.nas_score_sizes <- NULL
load_nas_score_sample_sizes <- function() {
  if (!is.null(.nas_score_sizes)) return(.nas_score_sizes)
  f <- file.path(SIGS, "nas_score_sample_sizes.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .nas_score_sizes <<- dt
  dt
}

# --- Per-fibrosis-stage dream results (Script 14b) ---
.fib_stage_dream <- NULL
load_fibrosis_stage_dream <- function() {
  if (!is.null(.fib_stage_dream)) return(.fib_stage_dream)
  f <- file.path(SIGS, "fibrosis_stage_dream.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("adj.P.Val" %in% names(dt) && !"padj" %in% names(dt))
    setnames(dt, "adj.P.Val", "padj")
  dt <- add_symbols(dt, "gene")
  .fib_stage_dream <<- dt
  dt
}

# --- Per-fibrosis-stage sample sizes ---
.fib_stage_sizes <- NULL
load_fibrosis_stage_sample_sizes <- function() {
  if (!is.null(.fib_stage_sizes)) return(.fib_stage_sizes)
  f <- file.path(SIGS, "fibrosis_stage_sample_sizes.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .fib_stage_sizes <<- dt
  dt
}

# --- Clean per-sample staged metadata (guards coarse/dropped cohorts) ------------
# Use this INSTEAD of reading unified_metadata.csv + binning fibrosis_stage directly.
# Returns per-sample metadata for the true-Kleiner-staged set with a ready `stage`
# factor (F0..F4), excluding COARSE_STAGE_COHORTS (Chen) + DROPPED_COHORTS (Gerhard).
.unified_meta_cache <- NULL
load_unified_metadata <- function() {
  if (!is.null(.unified_meta_cache)) return(.unified_meta_cache)
  f <- file.path(INT_META, "unified_metadata.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  .unified_meta_cache <<- fread(f)
  .unified_meta_cache
}

staged_disease_meta <- function(meta = NULL, disease_only = TRUE) {
  if (is.null(meta)) meta <- load_unified_metadata()
  m <- copy(as.data.table(meta))
  m[, .fib := suppressWarnings(as.integer(fibrosis_stage))]
  m <- m[!dataset %in% c(COARSE_STAGE_COHORTS, DROPPED_COHORTS)]
  m <- m[!is.na(.fib) & .fib %in% 0:4]
  if (disease_only && "group_binary" %in% names(m)) m <- m[group_binary != "Control"]
  m[, stage := factor(paste0("F", .fib), levels = paste0("F", 0:4))]
  m[, .fib := NULL]
  m[]
}

# NAS-stratified per-sample metadata. Drops NA-NAS (Bril's 26 unscored patients +
# Chen's wholly-absent NAS) and dropped cohorts; `nas_num` is the numeric score.
clean_nas_meta <- function(meta = NULL) {
  if (is.null(meta)) meta <- load_unified_metadata()
  m <- copy(as.data.table(meta))
  m[, nas_num := suppressWarnings(as.numeric(nas_score))]
  m <- m[!dataset %in% DROPPED_COHORTS][!is.na(nas_num)]
  m[]
}

# --- Per-fibrosis-stage GSEA results (Script 14c) ---
.fib_stage_gsea <- NULL
load_fibrosis_stage_gsea <- function() {
  if (!is.null(.fib_stage_gsea)) return(.fib_stage_gsea)
  f <- file.path(SIGS, "fibrosis_stage_gsea.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .fib_stage_gsea <<- dt
  dt
}

# --- NAS consecutive pairwise dream results (Script 14d) ---
.nas_consec_cache <- NULL
load_nas_consecutive <- function() {
  if (!is.null(.nas_consec_cache)) return(.nas_consec_cache)
  f <- file.path(SIGS, "nas_consecutive_dream.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("gene" %in% names(dt)) dt <- add_symbols(dt, "gene")
  .nas_consec_cache <<- dt
  dt
}

# --- Fibrosis consecutive pairwise dream results (Script 14d) ---
.fib_consec_cache <- NULL
load_fibrosis_consecutive <- function() {
  if (!is.null(.fib_consec_cache)) return(.fib_consec_cache)
  f <- file.path(SIGS, "fibrosis_consecutive_dream.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("gene" %in% names(dt)) dt <- add_symbols(dt, "gene")
  .fib_consec_cache <<- dt
  dt
}

# --- Progression landscape cells (Script 14d) ---
.landscape_cells_cache <- NULL
load_landscape_cells <- function() {
  if (!is.null(.landscape_cells_cache)) return(.landscape_cells_cache)
  f <- file.path(SIGS, "progression_landscape_cells.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .landscape_cells_cache <<- dt
  dt
}

# --- Progression landscape samples (Script 14d) ---
.landscape_samples_cache <- NULL
load_landscape_samples <- function() {
  if (!is.null(.landscape_samples_cache)) return(.landscape_samples_cache)
  f <- file.path(SIGS, "progression_landscape_samples.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .landscape_samples_cache <<- dt
  dt
}

# --- Sex-specific fgsea (male) ---
.sex_fgsea_male_cache <- NULL
load_sex_fgsea_male <- function() {
  if (!is.null(.sex_fgsea_male_cache)) return(.sex_fgsea_male_cache)
  f <- file.path(INT_RESULTS, "sex_fgsea_male.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .sex_fgsea_male_cache <<- dt
  dt
}

# --- Sex-specific fgsea (female) ---
.sex_fgsea_female_cache <- NULL
load_sex_fgsea_female <- function() {
  if (!is.null(.sex_fgsea_female_cache)) return(.sex_fgsea_female_cache)
  f <- file.path(INT_RESULTS, "sex_fgsea_female.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .sex_fgsea_female_cache <<- dt
  dt
}

# --- Per-sample dysregulation scores (Script 14e) ---
.dysreg_cache <- NULL
load_per_sample_dysregulation <- function() {
  if (!is.null(.dysreg_cache)) return(.dysreg_cache)
  f <- file.path(SIGS, "per_sample_dysregulation.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .dysreg_cache <<- dt
  dt
}

# --- NAS score GSEA (Script 14f) ---
.nas_gsea_cache <- NULL
load_nas_score_gsea <- function() {
  if (!is.null(.nas_gsea_cache)) return(.nas_gsea_cache)
  f <- file.path(SIGS, "nas_score_gsea.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .nas_gsea_cache <<- dt
  dt
}

# --- NAS x Fibrosis sample distribution grid (Script 14d) ---
.sample_dist_cache <- NULL
load_sample_distribution <- function() {
  if (!is.null(.sample_dist_cache)) return(.sample_dist_cache)
  f <- file.path(SIGS, "progression_sample_distribution.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .sample_dist_cache <<- dt
  dt
}

# --- Sex x disease interaction dream results ---
.sex_int_cache <- NULL
load_sex_interaction <- function() {
  if (!is.null(.sex_int_cache)) return(.sex_int_cache)
  f <- file.path(INT_RESULTS, "sex_interaction_dream.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  if ("gene" %in% names(dt)) dt <- add_symbols(dt, "gene")
  .sex_int_cache <<- dt
  dt
}

# --- Protein-transcript correlation ---
.prot_corr_cache <- NULL
load_protein_correlation <- function() {
  if (!is.null(.prot_corr_cache)) return(.prot_corr_cache)
  f <- file.path(PROTEOMICS_DIR, "protein_transcript_correlation.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .prot_corr_cache <<- dt
  dt
}

# --- Protein-validated drug targets ---
.prot_targets_cache <- NULL
load_protein_targets <- function() {
  if (!is.null(.prot_targets_cache)) return(.prot_targets_cache)
  f <- file.path(PROTEOMICS_DIR, "protein_validated_targets.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .prot_targets_cache <<- dt
  dt
}

# --- SCENIC+ disease regulons (24 TFs) ---
.regulon_cache <- NULL
load_disease_regulons <- function() {
  if (!is.null(.regulon_cache)) return(.regulon_cache)
  f <- file.path(ATAC_DIR, "scenic_plus", "disease_regulons.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .regulon_cache <<- dt
  dt
}

# --- SCENIC+ hepatocyte regulon targets (139 TF-target links) ---
.hep_regulon_cache <- NULL
load_hepatocyte_regulons <- function() {
  if (!is.null(.hep_regulon_cache)) return(.hep_regulon_cache)
  f <- file.path(ATAC_DIR, "scenic_plus", "hepatocyte_regulons.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .hep_regulon_cache <<- dt
  dt
}

# --- chromVAR TF motif activity (hepatocyte subset) ---
.chromvar_cache <- NULL
load_chromvar_hepatocyte <- function() {
  if (!is.null(.chromvar_cache)) return(.chromvar_cache)
  # A6 pseudoreplication fix (2026-06-20): serve the DONOR-LEVEL limma table
  # (chromvar_limma_per_ct.csv: cell_type/TF/logFC/adj.P.Val; 110 sig genome-wide,
  # 0 in hepatocytes) NOT the per-CELL Mann-Whitney table (chromvar_tf_activity.csv;
  # 4,832 "sig" = n-of-cells inflation). fig3_epigenomic_panels.R Panel(h) has a
  # schema shim that renames TF/logFC/adj.P.Val and zero-fills mean_deviation_*.
  f <- file.path(ATAC_DIR, "results", "chromvar_v2", "chromvar_limma_per_ct.csv")
  if (!file.exists(f))
    f <- file.path(ATAC_DIR, "results", "chromvar_v2", "chromvar_tf_activity.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  # Filter to hepatocyte rows (handle both naming conventions)
  if ("cell_type" %in% names(dt))
    dt <- dt[cell_type %in% c("Hepatocyte", "Hepatocytes")]
  .chromvar_cache <<- dt
  dt
}

# --- Spatial SVGs (condition comparison) ---
.svg_cache <- NULL
load_spatial_svgs <- function() {
  if (!is.null(.svg_cache)) return(.svg_cache)
  f <- file.path(SPATIAL_DIR, "svg", "differential_svgs.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .svg_cache <<- dt
  dt
}

# --- Spatial enrichment tests ---
.spatial_enr_cache <- NULL
load_spatial_enrichment <- function() {
  if (!is.null(.spatial_enr_cache)) return(.spatial_enr_cache)
  f <- file.path(SPATIAL_DIR, "integration", "spatial_enrichment_tests.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .spatial_enr_cache <<- dt
  dt
}

# --- GWAS-spatial convergence ---
.gwas_spatial_cache <- NULL
load_spatial_gwas <- function() {
  if (!is.null(.gwas_spatial_cache)) return(.gwas_spatial_cache)
  f <- file.path(SPATIAL_DIR, "gwas_spatial", "gwas_spatial_enrichment.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .gwas_spatial_cache <<- dt
  dt
}

# --- Drug validation table (10 drugs x evidence) ---
.drug_val_cache <- NULL
load_drug_validation_table <- function() {
  if (!is.null(.drug_val_cache)) return(.drug_val_cache)
  f <- file.path(DRUG, "clinical_drug_validation_table.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .drug_val_cache <<- dt
  dt
}

# --- Protein differential results (v2 corrected contrasts) ---
.prot_diff_cache <- NULL
load_protein_differential <- function() {
  if (!is.null(.prot_diff_cache)) return(.prot_diff_cache)
  # Prefer v3 (multi-contrast); fall back to v2, then v1
  f <- file.path(PROTEOMICS_DIR, "protein_differential_results_v3.csv")
  if (!file.exists(f)) f <- file.path(PROTEOMICS_DIR, "protein_differential_results_v2.csv")
  if (!file.exists(f)) f <- file.path(PROTEOMICS_DIR, "protein_differential_results.csv")
  if (!file.exists(f)) { message("WARNING: protein differential results not found"); return(NULL) }
  dt <- fread(f)
  .prot_diff_cache <<- dt
  dt
}

# --- Protein-transcript concordance (v3 multi-contrast; fall back to v2) ---
.prot_conc_v2_cache <- NULL
load_protein_concordance_v2 <- function() {
  if (!is.null(.prot_conc_v2_cache)) return(.prot_conc_v2_cache)
  f <- file.path(PROTEOMICS_DIR, "protein_transcript_concordance_v3.csv")
  if (!file.exists(f)) f <- file.path(PROTEOMICS_DIR, "protein_transcript_concordance_v2.csv")
  if (!file.exists(f)) { message("WARNING: protein concordance not found"); return(NULL) }
  dt <- fread(f)
  .prot_conc_v2_cache <<- dt
  dt
}

# --- Protein ranked enrichment (fgsea) ---
.prot_fgsea_cache <- NULL
load_protein_ranked_enrichment <- function() {
  if (!is.null(.prot_fgsea_cache)) return(.prot_fgsea_cache)
  f <- file.path(PROTEOMICS_DIR, "protein_ranked_enrichment.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .prot_fgsea_cache <<- dt
  dt
}

# --- Protein effect-size stratified detection ---
.prot_effsize_cache <- NULL
load_protein_effectsize_detection <- function() {
  if (!is.null(.prot_effsize_cache)) return(.prot_effsize_cache)
  f <- file.path(PROTEOMICS_DIR, "protein_effectsize_detection.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .prot_effsize_cache <<- dt
  dt
}

# --- Protein validation summary ---
.prot_valsummary_cache <- NULL
load_protein_validation_summary <- function() {
  if (!is.null(.prot_valsummary_cache)) return(.prot_valsummary_cache)
  f <- file.path(PROTEOMICS_DIR, "protein_validation_summary.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .prot_valsummary_cache <<- dt
  dt
}

# --- Spatial niche targets ---
.niche_cache <- NULL
load_spatial_niche_targets <- function() {
  if (!is.null(.niche_cache)) return(.niche_cache)
  f <- file.path(SPATIAL_DIR, "drug_repurposing", "spatial_niche_targets.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .niche_cache <<- dt
  dt
}

# --- Zonation classification ---
.zonation_cache <- NULL
load_zonation_classification <- function() {
  if (!is.null(.zonation_cache)) return(.zonation_cache)
  f <- file.path(BASE, "RNA-seq/results/zonation/deg_zonation_classification.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .zonation_cache <<- dt
  dt
}

# --- Single-cell pseudobulk DE (per cell type) ---
SC_PSEUDOBULK_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")

load_pseudobulk_de <- function(cell_type = "Hepatocytes") {
  f <- file.path(SC_PSEUDOBULK_DIR, paste0(cell_type, "_de.csv"))
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  dt
}

load_pseudobulk_all <- function() {
  if (!dir.exists(SC_PSEUDOBULK_DIR)) {
    message("WARNING: ", SC_PSEUDOBULK_DIR, " not found")
    return(NULL)
  }
  files <- list.files(SC_PSEUDOBULK_DIR, pattern = "_de\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  dt <- rbindlist(lapply(files, fread), fill = TRUE)
  dt
}

# --- Single-cell pathway enrichment (per cell type) ---
load_celltype_pathway_enrichment <- function() {
  f <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data/celltype_pathway_enrichment.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  dt
}

# --- Sex subsampling stability ---
.sex_sub_cache <- NULL
load_sex_subsampling <- function() {
  if (!is.null(.sex_sub_cache)) return(.sex_sub_cache)
  f <- file.path(AUDIT, "sex_subsampling", "sex_subsampling_fraction_summary.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .sex_sub_cache <<- dt
  dt
}

.sex_sub_iter_cache <- NULL
load_sex_subsampling_iter <- function() {
  if (!is.null(.sex_sub_iter_cache)) return(.sex_sub_iter_cache)
  f <- file.path(AUDIT, "sex_subsampling", "sex_subsampling_iter_metrics.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .sex_sub_iter_cache <<- dt
  dt
}

.sex_sub_gene_cache <- NULL
load_sex_subsampling_genes <- function() {
  if (!is.null(.sex_sub_gene_cache)) return(.sex_sub_gene_cache)
  f <- file.path(AUDIT, "sex_subsampling", "sex_subsampling_per_gene.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .sex_sub_gene_cache <<- dt
  dt
}

.sex_perm_cache <- NULL
load_sex_permutation <- function() {
  if (!is.null(.sex_perm_cache)) return(.sex_perm_cache)
  f <- file.path(AUDIT, "sex_subsampling", "sex_permutation_summary.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  dt <- fread(f)
  .sex_perm_cache <<- dt
  dt
}

load_sex_class_switching <- function() {
  f <- file.path(AUDIT, "sex_subsampling", "sex_subsampling_class_switching.csv")
  if (!file.exists(f)) { message("WARNING: ", f, " not found"); return(NULL) }
  fread(f)
}

message("load_figure_data.R loaded: ", length(ls(pattern = "^load_")), " loaders available")
